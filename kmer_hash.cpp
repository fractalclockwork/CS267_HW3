#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>
#include <fstream>
#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "butil.hpp"

// Initialize k-mers in Distributed HashMap
void initialize_kmers(DistributedHashMap &hashmap, const std::vector<kmer_pair> &kmers, 
                      std::vector<kmer_pair> &start_nodes) {
    for (const auto &kmer : kmers) {
        hashmap.insert(kmer.kmer_str(), kmer);
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }
    upcxx::barrier();
}

// Assemble contigs from start nodes using Distributed HashMap
std::list<std::list<kmer_pair>> assemble_contigs(DistributedHashMap &hashmap, 
                                                 const std::vector<kmer_pair> &start_nodes) {
    std::list<std::list<kmer_pair>> contigs;

    for (const auto &start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);

        while (contig.back().forwardExt() != 'F') {
            kmer_pair found;
            bool success = hashmap.find(contig.back().next_kmer().get(), found);
            if (!success) {
                throw std::runtime_error("Error: k-mer not found in Distributed HashMap.");
            }
            contig.push_back(found);
        }
        contigs.push_back(contig);
    }

    return contigs;
}

// Output results with metrics
void output_results(const std::list<std::list<kmer_pair>> &contigs, const std::string &test_prefix, 
                    int rank_id, double insert_time, double read_time, double total_time) {
    std::ofstream fout(test_prefix + "_" + std::to_string(rank_id) + ".dat");
    for (const auto &contig : contigs) {
        fout << extract_contig(contig) << std::endl;
    }

    BUtil::print("Rank %d reconstructed %d contigs with %d nodes."
                 " (%lf read, %lf insert, %lf total)\n",
                 rank_id, contigs.size(), std::accumulate(contigs.begin(), contigs.end(), 0, 
                 [](int sum, const std::list<kmer_pair> &contig) { return sum + contig.size(); }),
                 read_time, insert_time, total_time);
}

// Main function
int main(int argc, char **argv) {
    upcxx::init();

    int rank_id = upcxx::rank_me();
    int world_size = upcxx::rank_n();
    std::string kmer_fname = argv[1];
    std::string run_type = (argc >= 3) ? argv[2] : "";
    std::string test_prefix = (argc >= 4 && run_type == "test") ? argv[3] : "test";

    size_t hash_table_size = line_count(kmer_fname) * 2;
    DistributedHashMap hashmap(hash_table_size, rank_id, world_size);
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, world_size, rank_id);
    std::vector<kmer_pair> start_nodes;

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %lu for %lu kmers.\n", hash_table_size, kmers.size());
        BUtil::print("Finished reading kmers.\n");
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    initialize_kmers(hashmap, kmers, start_nodes);
    auto insert_time = std::chrono::high_resolution_clock::now();

    auto contigs = assemble_contigs(hashmap, start_nodes);
    auto end_time = std::chrono::high_resolution_clock::now();

    double insert_duration = std::chrono::duration<double>(insert_time - start_time).count();
    double read_duration = std::chrono::duration<double>(end_time - insert_time).count();
    double total_duration = std::chrono::duration<double>(end_time - start_time).count();

    if (run_type != "test") {
        BUtil::print("Finished inserting in %lf sec\n", insert_duration);
        BUtil::print("Assembled in %lf total\n", total_duration);
    }

    output_results(contigs, test_prefix, rank_id, insert_duration, read_duration, total_duration);
    upcxx::finalize();
    return 0;
}
