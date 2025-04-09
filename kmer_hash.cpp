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

#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "butil.hpp"

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
// Function: async_initialize_kmers
// Inserts k-mers asynchronously with dynamic batching.
//------------------------------------------------------
upcxx::future<> async_initialize_kmers(DistributedHashMap<std::string, kmer_pair> &hashmap,
                                        const std::vector<kmer_pair> &kmers,
                                        std::vector<kmer_pair> &start_nodes,
                                        size_t batch_size) {
    std::vector<upcxx::future<>> insert_futures;
    std::vector<std::pair<std::string, kmer_pair>> batch;

    for (size_t i = 0; i < kmers.size(); i++) {
        batch.emplace_back(kmers[i].kmer_str(), kmers[i]);
        if (kmers[i].backwardExt() == 'F') start_nodes.push_back(kmers[i]);

        if (batch.size() >= batch_size || i == kmers.size() - 1) {
            insert_futures.push_back(hashmap.batch_insert(batch));
            batch.clear();
        }
    }
    return upcxx::when_all(insert_futures).then([](auto){ return upcxx::make_future(); });
}

//------------------------------------------------------
// Function: assemble_contigs
// Traverses the hash map to reconstruct contigs.
//------------------------------------------------------
std::list<std::list<kmer_pair>> assemble_contigs(DistributedHashMap<std::string, kmer_pair> &hashmap, 
                                                 const std::vector<kmer_pair> &start_nodes, 
                                                 size_t &batch_size) {
    std::list<std::list<kmer_pair>> contigs;
    size_t success_count = 0, failure_count = 0;

    for (const auto &start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);
        while (contig.back().forwardExt() != 'F') {
            kmer_pair found;
            bool success = hashmap.find(contig.back().next_kmer().get(), found);
            if (!success) {
                failure_count++;
                break;
            }
            success_count++;
            contig.push_back(found);
        }
        contigs.push_back(contig);
    }
    
    batch_size = std::min(batch_size + size_t(50), size_t(1000));
    return contigs;
}
//------------------------------------------------------
// Function: output_results
// Writes each assembled contig to a file named based on test prefix and rank id.
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
    size_t n_kmers, hash_table_size, batch_size = 100;

    initialize_upcxx(argc, argv, kmer_fname, run_type, test_prefix, 
                     ks, n_kmers, hash_table_size, rank_id, world_size);

    DistributedHashMap<std::string, kmer_pair> hashmap(hash_table_size, rank_id, world_size);
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, world_size, rank_id);
    std::vector<kmer_pair> start_nodes;

    BUtil::print("Initializing hash table of size %lu for %lu kmers.\n", hash_table_size, n_kmers);

    auto start_time = std::chrono::high_resolution_clock::now();
    auto ins_fut = async_initialize_kmers(hashmap, kmers, start_nodes, batch_size);
    ins_fut.wait();
    auto end_insert = std::chrono::high_resolution_clock::now();

    auto contigs = assemble_contigs(hashmap, start_nodes, batch_size);
    auto end_read = std::chrono::high_resolution_clock::now();
    auto end_time = std::chrono::high_resolution_clock::now();

    upcxx::barrier();

    std::chrono::duration<double> insert = end_insert - start_time;
    std::chrono::duration<double> read = end_read - end_insert;
    std::chrono::duration<double> total = end_time - start_time;

    int numKmers = std::accumulate(contigs.begin(), contigs.end(), 0,
        [](int sum, const std::list<kmer_pair>& contig) { return sum + contig.size(); });

    BUtil::print("Finished inserting in %lf\n", insert.count());
    BUtil::print("Assembled in %lf total\n", total.count());

    printf("Rank %d reconstructed %lu contigs with %d nodes from %lu start nodes."
           " (%lf read, %lf insert, %lf total)\n",
           rank_id, contigs.size(), numKmers, start_nodes.size(), 
           read.count(), insert.count(), total.count());

    if (run_type == "test") {
        output_results(contigs, test_prefix, rank_id);
    }

    upcxx::finalize();
    return 0;
}
