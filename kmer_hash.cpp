#include <cstddef>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>
#include <chrono>
#include <fstream>
#include <stdexcept>
#include <string>

// Include necessary headers
#include "hash_map.hpp"   // Generic DistributedHashMap
#include "kmer_t.hpp"     // K-mer data structure
#include "read_kmers.hpp" // Handles k-mer file reading
#include "butil.hpp"      // Utility functions

//------------------------------------------------------
// Function: initialize_upcxx
// Parses command-line arguments, determines k-mer parameters, and sets up UPC++ values.
//------------------------------------------------------
void initialize_upcxx(int argc, char **argv, 
                      std::string &kmer_fname, std::string &run_type, std::string &test_prefix,
                      int &ks, size_t &n_kmers, size_t &hash_table_size, int &rank_id, int &world_size) {
    if (argc < 2) {
        BUtil::print("Usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }

    kmer_fname = argv[1];
    run_type = (argc >= 3) ? argv[2] : "";
    test_prefix = (run_type == "test" && argc >= 4) ? argv[3] : "test";

    ks = kmer_size(kmer_fname);
    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
                                 "-mers, while this binary is compiled for " +
                                 std::to_string(KMER_LEN) + "-mers.");
    }

    n_kmers = line_count(kmer_fname);
    hash_table_size = n_kmers * 2;
    rank_id = upcxx::rank_me();
    world_size = upcxx::rank_n();
}

//------------------------------------------------------
// Function: initialize_kmers
// Inserts k-mers into the distributed hash map and records "start nodes."
//------------------------------------------------------
/*
void initialize_kmers(DistributedHashMap<std::string, kmer_pair> &hashmap, 
                      const std::vector<kmer_pair> &kmers, 
                      std::vector<kmer_pair> &start_nodes) {
    for (const auto &kmer : kmers) {
        hashmap.insert(kmer.kmer_str(), kmer);
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }
    hashmap.process_requests();
}
*/
void initialize_kmers(DistributedHashMap<std::string, kmer_pair> &hashmap, 
    const std::vector<kmer_pair> &kmers, std::vector<kmer_pair> &start_nodes) {
    std::vector<std::pair<std::string, kmer_pair>> batched_kmers; // Collect k-mers in a batch

    for (const auto &kmer : kmers) {
        batched_kmers.emplace_back(kmer.kmer_str(), kmer);
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }

    // Use batch_insert to insert all k-mers efficiently
    hashmap.batch_insert(batched_kmers).wait();
    hashmap.process_requests();
}

//------------------------------------------------------
// Function: assemble_contigs
// Traverses the hash map from the start nodes, linking consecutive k-mers.
//------------------------------------------------------
std::list<std::list<kmer_pair>> assemble_contigs(DistributedHashMap<std::string, kmer_pair> &hashmap, 
                                                 const std::vector<kmer_pair> &start_nodes) {
    std::list<std::list<kmer_pair>> contigs;

    for (const auto &start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);

        while (contig.back().forwardExt() != 'F') {
            kmer_pair found;
            bool success = hashmap.find(contig.back().next_kmer().get(), found);
            if (!success) {
                throw std::runtime_error("Error: k-mer not found in hashmap.");
            }
            contig.push_back(found);
        }
        contigs.push_back(contig);
    }

    return contigs;
}

//------------------------------------------------------
// Function: output_results
// Writes each assembled contig to a file.
//------------------------------------------------------
void output_results(const std::list<std::list<kmer_pair>> &contigs, 
                    const std::string &test_prefix, int rank_id) {
    std::ofstream fout(test_prefix + "_" + std::to_string(rank_id) + ".dat");
    for (const auto &contig : contigs) {
        fout << extract_contig(contig) << std::endl;
    }
}

//------------------------------------------------------
// Main: Executes the k-mer processing pipeline.
//------------------------------------------------------
int main(int argc, char **argv) {
    upcxx::init();

    std::string kmer_fname, run_type, test_prefix;
    int ks, rank_id, world_size;
    size_t n_kmers, hash_table_size;

    initialize_upcxx(argc, argv, kmer_fname, run_type, test_prefix, 
                     ks, n_kmers, hash_table_size, rank_id, world_size);

    DistributedHashMap<std::string, kmer_pair> hashmap(hash_table_size, rank_id, world_size);
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, world_size, rank_id);
    std::vector<kmer_pair> start_nodes;

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %lu for %lu kmers.\n", hash_table_size, n_kmers);
        BUtil::print("Finished reading kmers.\n");
    }

    upcxx::barrier();
    auto start_time = std::chrono::high_resolution_clock::now();

    initialize_kmers(hashmap, kmers, start_nodes);
    upcxx::barrier();
    auto insert_time = std::chrono::high_resolution_clock::now();

    auto contigs = assemble_contigs(hashmap, start_nodes);
    upcxx::barrier();
    auto end_time = std::chrono::high_resolution_clock::now();

    // Reporting results
    if (run_type == "test") {
        output_results(contigs, test_prefix, rank_id);
    } else {
        BUtil::print("Finished inserting in %lf sec\n", 
                     std::chrono::duration<double>(insert_time - start_time).count());
        BUtil::print("Assembled in %lf sec\n", 
                     std::chrono::duration<double>(end_time - start_time).count());
    }

    upcxx::finalize();
    return 0;
}
