#include <cstddef>
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
        BUtil::print("Usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }

    std::string kmer_fname = argv[1];
    std::string run_type = (argc >= 3) ? argv[2] : "";
    std::string test_prefix = (run_type == "test" && argc >= 4) ? argv[3] : "test";

    int ks = kmer_size(kmer_fname);
    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
                                 "-mers, while this binary is compiled for " +
                                 std::to_string(KMER_LEN) + "-mers.");
    }

    size_t n_kmers = line_count(kmer_fname);
    size_t hash_table_size = n_kmers * (1.0 / 0.5);
    int rank_id = upcxx::rank_me();
    int world_size = upcxx::rank_n();
    
    DistributedHashMap hashmap(hash_table_size, rank_id, world_size);
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, world_size, rank_id);

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %lu for %lu kmers.\n", hash_table_size, n_kmers);
        BUtil::print("Finished reading kmers.\n");
    }

    upcxx::barrier();
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> start_nodes;
    for (auto& kmer : kmers) {
        hashmap.insert(kmer.kmer_str(), kmer).wait();
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }

    auto end_insert = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    
    auto start_read = std::chrono::high_resolution_clock::now();
    std::list<std::list<kmer_pair>> contigs;

    for (const auto& start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);
        while (contig.back().forwardExt() != 'F') {
            kmer_pair found;
            bool success = hashmap.find(contig.back().next_kmer().get(), found).wait();
            if (!success) {
                throw std::runtime_error("Error: k-mer not found in hashmap.");
            }
            contig.push_back(found);
        }
        contigs.push_back(contig);
    }

    auto end_read = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> insert_time = end_insert - start;
    std::chrono::duration<double> read_time = end_read - start_read;
    std::chrono::duration<double> total_time = end - start;

    int num_kmers = std::accumulate(contigs.begin(), contigs.end(), 0,
        [](int sum, const std::list<kmer_pair>& contig) { return sum + contig.size(); });

    if (run_type != "test") {
        BUtil::print("Finished inserting in %lf sec\n", insert_time.count());
        BUtil::print("Assembled in %lf sec total\n", total_time.count());
    }

    if (run_type == "verbose") {
        printf("Rank %d reconstructed %lu contigs with %d kmers from %lu start nodes."
               " (%lf sec read, %lf sec insert, %lf sec total)\n",
               rank_id, contigs.size(), num_kmers, start_nodes.size(),
               read_time.count(), insert_time.count(), total_time.count());
    }

    // **File Writing**
    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(rank_id) + ".dat");
        for (const auto& contig : contigs) {
            fout << extract_contig(contig) << std::endl;
        }
        fout.close();
    }

    upcxx::finalize();
    return 0;
}