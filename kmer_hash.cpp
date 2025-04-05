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
        BUtil::print("Usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        return 1;
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

    std::cout << "Rank " << rank_id << " processing " << kmers.size() << " kmers." << std::endl;

    auto start = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> start_nodes;
    for (auto& kmer : kmers) {
        hashmap.insert(kmer.kmer_str(), kmer).wait();

        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }

    upcxx::barrier();
    std::cout << "Rank " << rank_id << " identified " << start_nodes.size() << " start nodes." << std::endl;

    auto start_read = std::chrono::high_resolution_clock::now();
    std::list<std::list<kmer_pair>> contigs;

    for (const auto& start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);

        while (contig.back().forwardExt() != 'F') {
            kmer_pair found;
            bool success = hashmap.find(contig.back().next_kmer().get(), found).wait();

            if (!success) {
                std::cout << "Rank " << rank_id << " contig assembly failed at key: " << contig.back().kmer_str() << std::endl;
                break;
            }

            std::cout << "Rank " << rank_id << " extending contig with key: " << found.kmer_str() << std::endl;
            contig.push_back(found);
        }
        contigs.push_back(contig);
    }

    upcxx::barrier();
    auto end_read = std::chrono::high_resolution_clock::now();

    std::cout << "Rank " << rank_id << " assembled " << contigs.size() << " contigs." << std::endl;

    // Write contigs to file in append mode to prevent overwrites
    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(rank_id) + ".dat", std::ios::app);
        if (fout.fail()) {
            std::cerr << "Rank " << rank_id << ": Failed to open output file." << std::endl;
        } else {
            std::cout << "Rank " << rank_id << " writing " << contigs.size() << " contigs to file." << std::endl;
            for (const auto& contig : contigs) {
                fout << extract_contig(contig) << std::endl;
            }
            fout.flush();
            fout.close();
        }
    }

    upcxx::finalize();
    return 0;
}
