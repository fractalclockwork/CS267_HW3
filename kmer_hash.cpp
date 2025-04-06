#include <cstddef>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>
#include <chrono>
#include <fstream>

#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "butil.hpp"

// Function to initialize kmers into the Distributed HashMap
void initialize_kmers(DistributedHashMap &hashmap, 
    const std::vector<kmer_pair> &kmers, 
    std::vector<kmer_pair> &start_nodes) {
for (const auto &kmer : kmers) {
hashmap.insert(kmer.kmer_str(), kmer); // Insert each kmer one by one
if (kmer.backwardExt() == 'F') {
start_nodes.push_back(kmer);
}
}
}

// Function to assemble contigs from the start nodes using the Distributed HashMap
std::list<std::list<kmer_pair>> assemble_contigs(DistributedHashMap &hashmap, 
                               const std::vector<kmer_pair> &start_nodes) {
std::list<std::list<kmer_pair>> contigs;

for (const auto &start_kmer : start_nodes) {
std::list<kmer_pair> contig;
contig.push_back(start_kmer);

while (contig.back().forwardExt() != 'F') {
kmer_pair found;
bool success = hashmap.find(contig.back().next_kmer().get(), found);  // No .wait()
if (!success) {
throw std::runtime_error("Error: k-mer not found in hashmap.");
}
contig.push_back(found);
}
contigs.push_back(contig);
}

return contigs;
}

// Function to output the results after assembling the contigs
void output_results(const std::list<std::list<kmer_pair>> &contigs, 
                    const std::string &test_prefix, int rank_id) {
    std::ofstream fout(test_prefix + "_" + std::to_string(rank_id) + ".dat");
    for (const auto &contig : contigs) {
        fout << extract_contig(contig) << std::endl;
    }
}

// Function to initialize the UPC++ environment and handle command-line arguments
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

// Main function to execute the overall program
int main(int argc, char **argv) {
    upcxx::init();

    std::string kmer_fname, run_type, test_prefix;
    int ks, rank_id, world_size;
    size_t n_kmers, hash_table_size;

    // Initialize UPC++ environment and parse command-line arguments
    initialize_upcxx(argc, argv, kmer_fname, run_type, test_prefix, ks, n_kmers, hash_table_size, rank_id, world_size);

    // Initialize the Distributed HashMap
    DistributedHashMap hashmap(hash_table_size, rank_id, world_size);

    // Read the kmers from the file
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, world_size, rank_id);
    std::vector<kmer_pair> start_nodes;

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %lu for %lu kmers.\n", hash_table_size, n_kmers);
        BUtil::print("Finished reading kmers.\n");
    }

    upcxx::barrier();
    auto start = std::chrono::high_resolution_clock::now();

    // Initialize kmers into the Distributed HashMap
    initialize_kmers(hashmap, kmers, start_nodes);
    upcxx::barrier();
    auto insert_time = std::chrono::high_resolution_clock::now();

    // Assemble contigs from the start nodes
    auto contigs = assemble_contigs(hashmap, start_nodes);
    upcxx::barrier();
    auto end_time = std::chrono::high_resolution_clock::now();

    // Output results based on run type
    if (run_type == "test") {
        output_results(contigs, test_prefix, rank_id);
    }

    if (run_type != "test") {
        BUtil::print("Finished inserting in %lf sec\n", std::chrono::duration<double>(insert_time - start).count());
        BUtil::print("Assembled in %lf sec\n", std::chrono::duration<double>(end_time - start).count());
    }

    upcxx::finalize();
    return 0;
}
