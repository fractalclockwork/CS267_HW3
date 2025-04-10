#include <upcxx/upcxx.hpp>
#include <vector>
#include <list>
#include <chrono>
#include <fstream>
#include <stdexcept>
#include <string>
#include <numeric>
#include <cstdio>

#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "butil.hpp"

//------------------------------------------------------------------------
// initialize_kmer_processing:
// Parses command-line arguments and sets parameters for k-mer processing.
// This function preserves the original starter logic.
//------------------------------------------------------------------------
void initialize_kmer_processing(int argc, char **argv,
                                std::string &kmer_fname,
                                std::string &run_type,
                                std::string &test_prefix,
                                size_t &n_kmers,
                                size_t &hash_table_size) {
    if (argc < 2) {
        BUtil::print("usage: ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        exit(1);
    }
    kmer_fname = argv[1];
    run_type = (argc >= 3) ? argv[2] : "";
    test_prefix = (run_type == "test" && argc >= 4) ? argv[3] : "test";
    n_kmers = line_count(kmer_fname);
    hash_table_size = n_kmers * 2;  // For a load factor of 0.5.
}

//------------------------------------------------------------------------
// process_kmers:
// Inserts all k-mers into the hash map and records start nodes.
// Returns the insertion time.
//------------------------------------------------------------------------
double process_kmers(SerialHashMap &hashmap,
                     const std::vector<kmer_pair> &kmers,
                     std::vector<kmer_pair> &start_nodes) {
    auto start = std::chrono::high_resolution_clock::now();
    for (const auto &kmer : kmers) {
        if (!hashmap.insert(kmer))
            throw std::runtime_error("Error: HashMap is full!");
        if (kmer.backwardExt() == 'F')
            start_nodes.push_back(kmer);
    }
    auto end = std::chrono::high_resolution_clock::now();
    return std::chrono::duration<double>(end - start).count();
}

//------------------------------------------------------------------------
// assemble_contigs:
// Assembles contigs by traversing the hash map from each start node.
//------------------------------------------------------------------------
std::list<std::list<kmer_pair>> assemble_contigs(SerialHashMap &hashmap,
                                                 const std::vector<kmer_pair> &start_nodes) {
    std::list<std::list<kmer_pair>> contigs;
    for (const auto &start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);
        while (contig.back().forwardExt() != 'F') {
            kmer_pair next;
            if (!hashmap.find(contig.back().next_kmer(), next))
                throw std::runtime_error("Error: k-mer not found in hashmap.");
            contig.push_back(next);
        }
        contigs.push_back(contig);
    }
    return contigs;
}

//------------------------------------------------------------------------
// output_results:
// Writes each contig to an output file whose name includes the rank.
//------------------------------------------------------------------------
void output_results(const std::list<std::list<kmer_pair>> &contigs,
                    const std::string &filename_prefix) {
    // For serial execution using UPC++ we still use the rank (will be 0).
    std::ofstream fout(filename_prefix + "_" + std::to_string(upcxx::rank_me()) + ".dat");
    for (const auto &contig : contigs)
        fout << extract_contig(contig) << std::endl;
    fout.close();
}

int main(int argc, char **argv) {
    upcxx::init();  // Initialize UPC++ (this also prepares for future DHM extensions)

    std::string kmer_fname, run_type, test_prefix;
    size_t n_kmers, hash_table_size;
    initialize_kmer_processing(argc, argv, kmer_fname, run_type, test_prefix, n_kmers, hash_table_size);

    int ks = kmer_size(kmer_fname);
    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
                                   "-mers, while this binary is compiled for " +
                                   std::to_string(KMER_LEN) +
                                   "-mers. Modify packing.hpp and recompile.");
    }

    // Create the serial hash map.
    SerialHashMap hashmap(hash_table_size);
    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %zu for %zu kmers.\n", hash_table_size, n_kmers);
    }

    // Read and process k-mers.
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname);
    std::vector<kmer_pair> start_nodes;
    double insert_time = process_kmers(hashmap, kmers, start_nodes);
    BUtil::print("Finished inserting in %lf seconds.\n", insert_time);

    // Assemble contigs.
    auto start_read = std::chrono::high_resolution_clock::now();
    std::list<std::list<kmer_pair>> contigs = assemble_contigs(hashmap, start_nodes);
    auto end_read = std::chrono::high_resolution_clock::now();
    double read_time = std::chrono::duration<double>(end_read - start_read).count();

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
        output_results(contigs, test_prefix);
    }

    upcxx::finalize();
    return 0;
}
