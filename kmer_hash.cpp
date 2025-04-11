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
    size_t hash_table_size = static_cast<size_t>(n_kmers * (1.0 / 0.5));
    
    // Create the distributed hash table.
    HashMap hashmap(hash_table_size);

    // Read in the k-mers (each rank gets its portion).
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());
    upcxx::barrier(); // Synchronize before proceeding.

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> start_nodes;
    // Insert each k-mer into the distributed hash table.
    for (auto& kmer : kmers) {
        bool success = hashmap.insert(kmer);
        if (!success) {
            throw std::runtime_error("Error: HashMap insertion failed for k-mer " + kmer.kmer_str());
        }
        // Debug: if verbose mode is enabled, print a message for each inserted starting node.
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
            BUtil::print("Rank %d: Marking k-mer %s as a start node.\n",
                         upcxx::rank_me(), kmer.kmer_str().c_str());
        }
    }

    auto end_insert = std::chrono::high_resolution_clock::now();
    upcxx::barrier();

    double insert_time = std::chrono::duration<double>(end_insert - start).count();
    // Always print insertion metrics.
    BUtil::print("Finished inserting %zu k-mers in %lf seconds (Rank %d).\n",
                 kmers.size(), insert_time, upcxx::rank_me());
    upcxx::barrier();

    auto start_read = std::chrono::high_resolution_clock::now();

    // Assemble contigs by traversing from each starting node.
    std::list<std::list<kmer_pair>> contigs;
    for (const auto& start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);
        // Continue traversal until the forward extension is 'F'
        while (contig.back().forwardExt() != 'F') {
            kmer_pair next_kmer;
            bool success = hashmap.find(contig.back().next_kmer(), next_kmer);
            if (!success) {
                BUtil::print("Rank %d: Contig lookup failed at k-mer %s; failing contig assembly.\n",
                             upcxx::rank_me(), contig.back().next_kmer().get().c_str());
                throw std::runtime_error("Error: k-mer lookup failed in contig assembly.");
            }
            contig.push_back(next_kmer);
        }
        contigs.push_back(contig);
    }

    auto end_read = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> read_time = end_read - start_read;
    std::chrono::duration<double> total_time = end - start;

    // Print verbose reconstruction metrics.
    if (run_type == "verbose") {
        // Compute total nodes across all contigs.
        int numKmers = 0;
        for (const auto& contig : contigs)
            numKmers += contig.size();
        BUtil::print("Rank %d: Reconstructed %zu contigs with %d nodes from %zu start nodes. (read=%lf, insert=%lf, total=%lf seconds)\n",
                     upcxx::rank_me(), contigs.size(), numKmers, start_nodes.size(),
                     read_time.count(), insert_time, total_time.count());
    }
    
    // In test mode, write out the reconstructed contigs.
    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        if (!fout) {
            throw std::runtime_error("Error: Unable to open output file.");
        }
        for (const auto& contig : contigs) {
            fout << extract_contig(contig) << "\n";
        }
        fout.close();
    }

    upcxx::finalize();
    return 0;
}
