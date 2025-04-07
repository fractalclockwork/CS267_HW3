#pragma once
#include <unordered_map>
#include <string>
#include <vector>
#include <upcxx/upcxx.hpp>
#include "kmer_t.hpp"

class DistributedHashMap {
private:
    using dobj_map_t = upcxx::dist_object<std::unordered_map<std::string, kmer_pair>>;
    dobj_map_t local_map;

    // Local cache for frequently accessed k-mers
    std::unordered_map<std::string, kmer_pair> local_cache;

    size_t table_size_;
    int rank_id_;
    int world_size_;

    int get_target_rank(const std::string &key) const {
        return std::hash<std::string>{}(key) % world_size_;
    }

public:
    DistributedHashMap(size_t table_size, int rank_id, int world_size)
        : table_size_(table_size), rank_id_(rank_id), world_size_(world_size), local_map({}) {}

    // **Batch Insert Operation**
    upcxx::future<> batch_insert(const std::vector<std::pair<std::string, kmer_pair>> &entries) {
        std::unordered_map<int, std::vector<std::pair<std::string, kmer_pair>>> grouped_inserts;

        // **Step 1: Group entries by target rank**
        for (const auto &entry : entries) {
            grouped_inserts[get_target_rank(entry.first)].push_back(entry);
            local_cache[entry.first] = entry.second;  // Cache newly inserted k-mers
        }

        std::vector<upcxx::future<>> rpc_futures;

        // **Step 2: Perform insertions**
        for (const auto &[target_rank, batch] : grouped_inserts) {
            if (target_rank == rank_id_) {
                for (const auto &kv : batch) {
                    local_map->insert({kv.first, kv.second});
                }
            } else {
                rpc_futures.push_back(upcxx::rpc(target_rank,
                    [](dobj_map_t &lmap, const std::vector<std::pair<std::string, kmer_pair>> &batch) {
                        for (const auto &kv : batch) {
                            lmap->insert({kv.first, kv.second});
                        }
                    }, local_map, batch));
            }
        }

        if (rpc_futures.empty()) {
            return upcxx::make_future();
        }

        for (auto &future : rpc_futures) {
            future.wait();
        }

        return upcxx::make_future();
    }

    // **Updated Insert Function to Use Batch Insert**
    void insert(const std::string &key, const kmer_pair &value) {
        batch_insert({{key, value}}).wait();
    }

    // **Optimized Lookup Using Local Cache**
    bool find(const std::string &key, kmer_pair &result) {
        // **Step 1: Check Local Cache First**
        auto cache_it = local_cache.find(key);
        if (cache_it != local_cache.end()) {
            result = cache_it->second;
            return true;
        }

        int target_rank = get_target_rank(key);
        if (target_rank == rank_id_) {
            auto it = local_map->find(key);
            if (it == local_map->end()) {
                return false;
            }
            result = it->second;
            local_cache[key] = result;  // Cache retrieved value
            return true;
        } else {
            kmer_pair found_kmer = upcxx::rpc(target_rank,
                [](dobj_map_t &lmap, const std::string &key) -> kmer_pair {
                    auto it = lmap->find(key);
                    return (it != lmap->end()) ? it->second : kmer_pair();
                }, local_map, key).wait();

            if (!found_kmer.kmer_str().empty()) {
                result = found_kmer;
                local_cache[key] = result;  // Cache remote lookup result
                return true;
            }
            return false;
        }
    }

    void process_requests() {
        upcxx::progress(upcxx::progress_level::user);
    }
};