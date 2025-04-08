#pragma once
#include <unordered_map>
#include <vector>
#include <upcxx/upcxx.hpp>

//======================================================
// Class: DistributedHashMap
// Implements a **generic distributed hash table** that can store and retrieve 
// any key-value pair across UPC++ distributed memory.
//======================================================
template<typename Key, typename Value>
class DistributedHashMap {
private:
    using dist_map_t = upcxx::dist_object<std::unordered_map<Key, Value>>;
    dist_map_t local_map;

    // Local cache for frequently accessed data
    std::unordered_map<Key, Value> local_cache;

    size_t table_size_;  // Logical hash table size (maintained for interface consistency)
    int rank_id_;        // Current UPC++ rank identifier
    int world_size_;     // Total number of UPC++ ranks in execution

    //======================================================
    // Function: get_target_rank
    // Hash function to determine which rank stores a given key.
    //======================================================
    int get_target_rank(const Key &key) const {
        return std::hash<Key>{}(key) % world_size_;
    }

public:
    //======================================================
    // Constructor: DistributedHashMap
    // Initializes the distributed hash table.
    //======================================================
    DistributedHashMap(size_t table_size, int rank_id, int world_size)
        : table_size_(table_size), rank_id_(rank_id), world_size_(world_size), local_map({}) {}

    //======================================================
    // Function: batch_insert
    // Efficiently inserts multiple key-value pairs into the distributed hash table.
    //======================================================
    upcxx::future<> batch_insert(const std::vector<std::pair<Key, Value>> &entries) {
        std::unordered_map<int, std::vector<std::pair<Key, Value>>> grouped_inserts;

        for (const auto &entry : entries) {
            grouped_inserts[get_target_rank(entry.first)].push_back(entry);
            local_cache[entry.first] = entry.second;  // Cache inserted entries locally
        }

        std::vector<upcxx::future<>> rpc_futures;

        for (const auto &[target_rank, batch] : grouped_inserts) {
            if (target_rank == rank_id_) {
                for (const auto &kv : batch) {
                    local_map->insert({kv.first, kv.second});
                }
            } else {
                rpc_futures.push_back(upcxx::rpc(target_rank,
                    [](dist_map_t &lmap, const std::vector<std::pair<Key, Value>> &batch) {
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

    //======================================================
    // Function: insert
    // Wrapper for batch insertion of a single key-value pair.
    //======================================================
    void insert(const Key &key, const Value &value) {
        batch_insert({{key, value}}).wait();
    }

    //======================================================
    // Function: find
    // Retrieves a value from the distributed hash table using a key.
    // Returns true if the value was found, false otherwise.
    //======================================================
    bool find(const Key &key, Value &result) {
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
            local_cache[key] = result;  // Cache locally
            return true;
        } else {
            Value found_value = upcxx::rpc(target_rank,
                [](dist_map_t &lmap, const Key &key) -> Value {
                    auto it = lmap->find(key);
                    return (it != lmap->end()) ? it->second : Value{};
                }, local_map, key).wait();

            if (!(found_value == Value{})) {
                result = found_value;
                local_cache[key] = result;  // Cache remotely retrieved value
                return true;
            }
            return false;
        }
    }

    //======================================================
    // Function: process_requests
    // Executes pending UPC++ asynchronous operations.
    //======================================================
    void process_requests() {
        upcxx::progress(upcxx::progress_level::user);
    }
};
