#include <upcxx/upcxx.hpp>
#include <vector>
#include <list>
#include <chrono>
#include <fstream>
#include <stdexcept>
#include <string>
#include <numeric>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "butil.hpp"

//---------------------------------------------------------------------
// Helper Functions for Parallel de Bruijn Graph Traversal
//---------------------------------------------------------------------

// Convert a pkmer_t to a std::string via its get() method.
std::string pkmer_to_string(const pkmer_t &pk) {
    return pk.get();
}

// Compute the next k-mer for forward traversal.
// Given k-mer kp, drop its first character and append its forward extension.
pkmer_t compute_next_kmer(const kmer_pair &kp) {
    std::string curr = pkmer_to_string(kp.kmer);
    if (curr.empty())
        throw std::runtime_error("Empty k-mer in compute_next_kmer");
    // Next k-mer = curr.substr(1) + forward extension char.
    std::string next_str = curr.substr(1) + kp.forwardExt();
    return pkmer_t(next_str);
}

// Compute the previous k-mer for backward traversal.
// Given k-mer kp, drop its last character and prepend its backward extension.
pkmer_t compute_prev_kmer(const kmer_pair &kp) {
    std::string curr = pkmer_to_string(kp.kmer);
    if (curr.empty())
        throw std::runtime_error("Empty k-mer in compute_prev_kmer");
    // Previous k-mer = backward extension char + curr.substr(0, k-1)
    std::string prev_str = std::string(1, kp.backwardExt()) + curr.substr(0, curr.size()-1);
    return pkmer_t(prev_str);
}

// Extend contig forward starting from seed.
std::string extend_forward(const HashMap &hm, const kmer_pair &seed) {
    kmer_pair cur = seed;
    std::string contig = pkmer_to_string(cur.kmer);
    while (cur.forwardExt() != 'F') { // 'F' indicates no forward edge.
        char f = cur.forwardExt();
        pkmer_t next_key = compute_next_kmer(cur);
        kmer_pair next;
        bool found = hm.find(next_key, next);
        if (!found)
            break;
        contig.push_back(f);
        cur = next;
    }
    return contig;
}

// Extend contig backward starting from seed.
std::string extend_backward(const HashMap &hm, const kmer_pair &seed) {
    kmer_pair cur = seed;
    std::string contig = pkmer_to_string(cur.kmer);
    while (cur.backwardExt() != 'F') { // 'F' indicates no backward edge.
        char b = cur.backwardExt();
        pkmer_t prev_key = compute_prev_kmer(cur);
        kmer_pair prev;
        bool found = hm.find(prev_key, prev);
        if (!found)
            break;
        contig.insert(contig.begin(), b);
        cur = prev;
    }
    return contig;
}

// Traverse the contig by merging backward and forward extensions.
// Both left and right include the seed; remove duplicate seed region.
std::string traverse_contig(const HashMap &hm, const kmer_pair &seed) {
    std::string left = extend_backward(hm, seed);
    std::string right = extend_forward(hm, seed);
    std::string seed_str = pkmer_to_string(seed.kmer);
    if (right.size() < seed_str.size())
        return right;
    return left + right.substr(seed_str.size());
}

//---------------------------------------------------------------------
// Initialize processing parameters.
void initialize_kmer_processing(int argc, char** argv, std::string &kmer_fname,
                                std::string &run_type, std::string &test_prefix,
                                size_t &n_kmers, size_t &hash_table_size) {
    if (argc < 2) {
        BUtil::print("usage: ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        exit(1);
    }
    kmer_fname = argv[1];
    run_type = (argc >= 3) ? argv[2] : "";
    test_prefix = (run_type == "test" && argc >= 4) ? argv[3] : "test";
    n_kmers = line_count(kmer_fname);
    hash_table_size = n_kmers * 2;  // Load factor 0.5
}

int main(int argc, char** argv) {
    upcxx::init();
    srand((unsigned)time(NULL));

    std::string kmer_fname, run_type, test_prefix;
    size_t n_kmers, hash_table_size;
    initialize_kmer_processing(argc, argv, kmer_fname, run_type, test_prefix, n_kmers, hash_table_size);

    int ks = kmer_size(kmer_fname);
    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
             "-mers, while this binary is compiled for " + std::to_string(KMER_LEN) +
             "-mers. Modify packing.hpp and recompile.");
    }

    // Mark overall start time.
    auto start_total = std::chrono::high_resolution_clock::now();

    // Create the serial hash map.
    HashMap hashmap(hash_table_size);
    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %zu for %zu kmers.\n", hash_table_size, n_kmers);
    }

    // Read all k-mers.
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname);

    // Insertion Phase.
    auto start_insert = std::chrono::high_resolution_clock::now();
    for (const auto &kmer : kmers) {
        if (!hashmap.insert(kmer))
            throw std::runtime_error("Error: HashMap is full!");
    }
    auto end_insert = std::chrono::high_resolution_clock::now();
    double insert_time = std::chrono::duration<double>(end_insert - start_insert).count();
    BUtil::print("Finished inserting in %lf seconds.\n", insert_time);

    // Identify start nodes for traversal.
    std::vector<kmer_pair> start_nodes;
    for (const auto &kmer : kmers) {
        if (kmer.backwardExt() == 'F') {  // Has no valid incoming edge.
            start_nodes.push_back(kmer);
        }
    }
    if (start_nodes.empty()) {
        // If none found, fall back to a random seed.
        start_nodes.push_back(kmers[rand() % kmers.size()]);
    }

    // Traversal Phase: Assemble contigs for every start node.
    std::vector<std::string> contigs;
    auto start_traverse = std::chrono::high_resolution_clock::now();
    for (const auto &seed : start_nodes) {
        std::string contig = traverse_contig(hashmap, seed);
        contigs.push_back(contig);
    }
    auto end_traverse = std::chrono::high_resolution_clock::now();
    double trav_time = std::chrono::duration<double>(end_traverse - start_traverse).count();
    BUtil::print("Assembled %zu contigs in %lf seconds.\n", contigs.size(), trav_time);

    // Mark overall end time.
    auto end_total = std::chrono::high_resolution_clock::now();
    double total_time = std::chrono::duration<double>(end_total - start_total).count();
    BUtil::print("Total execution time: %lf seconds.\n", total_time);

    // In test mode, output all contigs to a file (one per line).
    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        for (const auto &c : contigs)
            fout << c << std::endl;
        fout.close();
    }

    upcxx::finalize();
    return 0;
}
