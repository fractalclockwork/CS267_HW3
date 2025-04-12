#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>
#include <fstream>
#include <sstream>
#include "hash_map.hpp"      // our refactored distributed hash table
#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "butil.hpp"

// -------------------------------------------------------------------------
// Function: initialize_kmers
//   Splits the local k-mers into a batch and inserts them in one call.
//   Also collects start nodes (k-mers with backward extension 'F').
void initialize_kmers(DistributedHashMap &hashmap, 
                      const std::vector<kmer_pair> &kmers, 
                      std::vector<kmer_pair> &start_nodes) {
    // Instead of per-insert RPC calls, we do a full batch insertion.
    hashmap.insert_all(kmers);
    // Also record start nodes.
    for (const auto &kmer : kmers) {
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }
    upcxx::barrier();
}

// -------------------------------------------------------------------------
// Function: assemble_contigs
//   Uses distributed find operations to follow forward extensions and build contigs.
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

// -------------------------------------------------------------------------
// Function: output_results
//   Prints metrics in the same format as the starter code and writes contigs.
void output_results(const std::list<std::list<kmer_pair>> &contigs, 
                    const std::string &test_prefix, int rank_id, 
                    double insert_time, double assembly_time, double total_time) {
    // Output assembled contigs (for test mode).
    std::ofstream fout(test_prefix + "_" + std::to_string(rank_id) + ".dat");
    for (const auto &contig : contigs) {
        fout << extract_contig(contig) << std::endl;
    }
    fout.close();
    
    // Print formatted metrics (as in the starter code).
    BUtil::print("Rank %d reconstructed %d contigs with %d nodes from %d start nodes. "
                 "(%lf read, %lf insert, %lf total)\n",
         rank_id,
         (int)contigs.size(),
         std::accumulate(contigs.begin(), contigs.end(), 0, 
            [](int sum, const std::list<kmer_pair>& contig) { return sum + contig.size(); }),
         0, // In this refactored design, start node count could be added here
         assembly_time, insert_time, total_time);
}

// -------------------------------------------------------------------------
// Main: Distributed genome assembler using a scalable distributed hash table.
// -------------------------------------------------------------------------
int main(int argc, char **argv) {
    upcxx::init();

    if(argc < 2){
        BUtil::print("Usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }
    
    std::string kmer_fname = std::string(argv[1]);
    std::string run_type = (argc >= 3) ? std::string(argv[2]) : "";
    std::string test_prefix = "test";
    if(run_type == "test" && argc >= 4){
        test_prefix = std::string(argv[3]);
    }
    
    int ks = kmer_size(kmer_fname);
    if(ks != KMER_LEN){
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
            "-mers, while this binary is compiled for " + std::to_string(KMER_LEN) +
            "-mers. Modify packing.hpp and recompile.");
    }
    
    size_t n_kmers = line_count(kmer_fname);
    // Load factor of 0.5 implies table size = n_kmers*2.
    size_t hash_table_size = n_kmers * 2;
    
    int rank_id = upcxx::rank_me();
    int world_size = upcxx::rank_n();
    
    if(run_type == "verbose"){
        BUtil::print("Initializing hash table of size %lu for %lu kmers.\n", hash_table_size, n_kmers);
    }
    
    // Create our scalable distributed hash map.
    DistributedHashMap hashmap(hash_table_size, rank_id, world_size);
    
    // Read the k-mers (each rank gets a portion).
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, world_size, rank_id);
    if(run_type == "verbose"){
        BUtil::print("Finished reading kmers.\n");
    }
    upcxx::barrier();
    
    // Timing: begin insertion.
    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<kmer_pair> start_nodes;
    initialize_kmers(hashmap, kmers, start_nodes);
    auto insert_time = std::chrono::high_resolution_clock::now();
    
    // Assemble contigs using distributed lookups.
    auto contigs = assemble_contigs(hashmap, start_nodes);
    upcxx::barrier();
    auto end_time = std::chrono::high_resolution_clock::now();
    
    double insert_duration = std::chrono::duration<double>(insert_time - start_time).count();
    double assembly_duration = std::chrono::duration<double>(end_time - insert_time).count();
    double total_duration = std::chrono::duration<double>(end_time - start_time).count();
    
    if(run_type != "test"){
        BUtil::print("Finished inserting in %lf sec\n", insert_duration);
        BUtil::print("Assembled in %lf total\n", total_duration);
    } else {
        output_results(contigs, test_prefix, rank_id, insert_duration, assembly_duration, total_duration);
    }
    
    upcxx::finalize();
    return 0;
}
