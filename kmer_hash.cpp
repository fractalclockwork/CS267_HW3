#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>

#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "butil.hpp"

int main(int argc, char** argv) {
    upcxx::init();

    if (argc < 2) {
        BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }

    std::string kmer_fname = std::string(argv[1]);
    std::string run_type = (argc >= 3) ? std::string(argv[2]) : "";

    std::string test_prefix = (run_type == "test" && argc >= 4) ? std::string(argv[3]) : "test";

    int ks = kmer_size(kmer_fname);
    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: Mismatched k-mer length. Modify packing.hpp and recompile.");
    }

    size_t n_kmers = line_count(kmer_fname);
    size_t hash_table_size = n_kmers * (1.0 / 0.5);
    
    HashMap hashmap(hash_table_size);

    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());
    upcxx::barrier();

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> start_nodes;
    for (auto& kmer : kmers) {
        bool success = hashmap.insert(kmer);
        if (!success) {
            throw std::runtime_error("Error: HashMap insertion failed!");
        }
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }

    auto end_insert = std::chrono::high_resolution_clock::now();
    upcxx::barrier();

    double insert_time = std::chrono::duration<double>(end_insert - start).count();
    BUtil::print("Finished inserting in %lf seconds\n", insert_time);
    upcxx::barrier();

    auto start_read = std::chrono::high_resolution_clock::now();

    std::list<std::list<kmer_pair>> contigs;
    for (const auto& start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);
        while (contig.back().forwardExt() != 'F') {
            kmer_pair next_kmer;
            bool success = hashmap.find(contig.back().next_kmer(), next_kmer);
            if (!success) {
                throw std::runtime_error("Error: k-mer lookup failed.");
            }
            contig.push_back(next_kmer);
        }
        contigs.push_back(contig);
    }

    auto end_read = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_time = end_read - start_read;
    std::chrono::duration<double> insert_duration = end_insert - start;
    std::chrono::duration<double> total_time = end - start;

    if (run_type != "test") {
        BUtil::print("Assembled in %lf seconds total\n", total_time.count());
    }

    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        for (const auto& contig : contigs) {
            fout << extract_contig(contig) << std::endl;
        }
        fout.close();
    }

    upcxx::finalize();
    return 0;
}
