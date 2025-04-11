#pragma once
#include <upcxx/upcxx.hpp>
#include <unordered_map>
#include <string>
#include "kmer_t.hpp"

// Distributed HashMap for parallel k-mer storage using UPC++ dist_object and RPC
class DistributedHashMap {
private:
    using dobj_map_t = upcxx::dist_object<std::unordered_map<std::string, kmer_pair>>;
    dobj_map_t local_map;
    size_t table_size_;
    int rank_id_;
    int world_size_;

    // Compute the target rank for a given key
    int get_target_rank(const std::string &key) const {
        return std::hash<std::string>{}(key) % world_size_;
    }

    // Insert locally within this rank's hash table
    void insert_locally(const std::string &key, const kmer_pair &value) {
        (*local_map)[key] = value;
    }

    // Perform remote insertion via UPC++ RPC
    void insert_remotely(int target_rank, const std::string &key, const kmer_pair &value) {
        upcxx::rpc(target_rank,
            [](upcxx::dist_object<std::unordered_map<std::string, kmer_pair>> &lmap, 
               const std::string &key, const kmer_pair &value) {
                (*lmap)[key] = value;
            },
            local_map, key, value).wait();
    }

    // Find locally within this rank's hash table
    bool find_locally(const std::string &key, kmer_pair &result) {
        auto it = local_map->find(key);
        if (it == local_map->end()) return false;
        result = it->second;
        return true;
    }

    // Perform remote find operation via UPC++ RPC
    bool find_remotely(int target_rank, const std::string &key, kmer_pair &result) {
        auto future_result = upcxx::rpc(target_rank,
            [](upcxx::dist_object<std::unordered_map<std::string, kmer_pair>> &lmap, 
               const std::string &key) -> upcxx::future<kmer_pair> {
                auto it = lmap->find(key);
                return upcxx::make_future((it != lmap->end()) ? it->second : kmer_pair());
            },
            local_map, key);
        
        kmer_pair found_kmer = future_result.wait();
        if (!found_kmer.kmer_str().empty()) {
            result = found_kmer;
            return true;
        }
        return false;
    }

public:
    // Constructor initializes the distributed hash map
    DistributedHashMap(size_t table_size, int rank_id, int world_size)
        : table_size_(table_size), rank_id_(rank_id), world_size_(world_size), local_map({}) {}

    // Insert a key-value pair into the map
    void insert(const std::string &key, const kmer_pair &value) {
        int target_rank = get_target_rank(key);
        if (target_rank == rank_id_) {
            insert_locally(key, value);
        } else {
            insert_remotely(target_rank, key, value);
        }
    }

    // Find a key in the map
    bool find(const std::string &key, kmer_pair &result) {
        int target_rank = get_target_rank(key);
        return (target_rank == rank_id_) ? find_locally(key, result) : find_remotely(target_rank, key, result);
    }

    // Synchronize processes and process incoming RPC calls
    void process_requests() {
        upcxx::progress(upcxx::progress_level::user);
        upcxx::barrier();
    }
};
