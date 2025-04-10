#include <upcxx/upcxx.hpp>
#include <vector>
#include <list>
#include <chrono>
#include <fstream>
#include <stdexcept>
#include <string>
#include <numeric>
#include <cstdio>
#include "hash_map.hpp"  // This file now contains DistributedHashMap
#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "butil.hpp"

int main(int argc, char **argv) {
    upcxx::init();

    if (argc < 2) {
        BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash_dhm kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }
    
    std::string kmer_fname(argv[1]);
    std::string run_type = (argc >= 3) ? argv[2] : "";
    std::string test_prefix = (run_type == "test" && argc >= 4) ? argv[3] : "test";
    
    int ks = kmer_size(kmer_fname);
    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
                                   "-mers, while this binary is compiled for " + std::to_string(KMER_LEN) +
                                   "-mers.");
    }
    
    size_t n_kmers = line_count(kmer_fname);
    size_t hash_table_size = n_kmers * 2;  // Maintain a 0.5 load factor.
    
    // Create the distributed hash map.
    DistributedHashMap<pkmer_t, kmer_pair> dhash(hash_table_size);
    if (run_type == "verbose") {
        BUtil::print("Initializing distributed hash table of size %zu for %zu kmers.\n", hash_table_size, n_kmers);
    }
    
    // Read k-mers (each rank reads its partition; adjust as needed).
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());
    std::vector<kmer_pair> start_nodes;
    auto insert_start = std::chrono::high_resolution_clock::now();
    
    // Insert each k-mer into the distributed hash map.
    for (const auto &kmer : kmers) {
        bool success = dhash.insert(kmer.kmer, kmer);
        if (!success) {
            throw std::runtime_error("Error: Distributed Hash Map is full!");
        }
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }
    
    auto insert_end = std::chrono::high_resolution_clock::now();
    double insert_time = std::chrono::duration<double>(insert_end - insert_start).count();
    BUtil::print("Finished inserting in %lf seconds.\n", insert_time);
    
    // Assemble contigs.
    std::list<std::list<kmer_pair>> contigs;
    auto read_start = std::chrono::high_resolution_clock::now();
    for (const auto &start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);
        while (contig.back().forwardExt() != 'F') {
            kmer_pair next;
            bool found = dhash.find(contig.back().next_kmer(), next);
            if (!found) {
                throw std::runtime_error("Error: k-mer not found in Distributed Hash Map.");
            }
            contig.push_back(next);
        }
        contigs.push_back(contig);
    }
    auto read_end = std::chrono::high_resolution_clock::now();
    double read_time = std::chrono::duration<double>(read_end - read_start).count();
    double total_time = insert_time + read_time;
    int numKmers = std::accumulate(
        contigs.begin(), contigs.end(), 0,
        [](int sum, const std::list<kmer_pair> &contig) {
            return sum + static_cast<int>(contig.size());
        }
    );
    
    BUtil::print("Assembled in %lf total seconds.\n", total_time);
    if (run_type == "verbose") {
        printf("Reconstructed %zu contigs with %d nodes from %zu start nodes. (%lf read, %lf insert, %lf total)\n",
               contigs.size(), numKmers, start_nodes.size(), read_time, insert_time, total_time);
    }
    
    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        for (const auto &contig : contigs)
            fout << extract_contig(contig) << std::endl;
        fout.close();
    }
    
    upcxx::finalize();
    return 0;
}
