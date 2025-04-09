#pragma once
#include <unordered_map>
#include <vector>
#include <upcxx/upcxx.hpp>

//======================================================
// Class: DistributedHashMap
// Implements a **generic distributed hash table** with dynamic batching and caching.
// This data structure is application agnostic.
//======================================================
template<typename Key, typename Value>
class DistributedHashMap {
private:
    using dist_map_t = upcxx::dist_object<std::unordered_map<Key, Value>>;
    dist_map_t local_map;
    std::unordered_map<Key, Value> local_cache;
    size_t table_size_;
    int rank_id_;
    int world_size_;

    //------------------------------------------------------
    // Function: get_target_rank
    // Determines which rank is responsible for the given key.
    //------------------------------------------------------
    int get_target_rank(const Key &key) const {
        return std::hash<Key>{}(key) % world_size_;
    }

public:
    //------------------------------------------------------
    // Constructor: DistributedHashMap
    // Initializes the distributed hash table.
    //------------------------------------------------------
    DistributedHashMap(size_t table_size, int rank_id, int world_size)
        : table_size_(table_size), rank_id_(rank_id), world_size_(world_size), local_map({}) {}

    //------------------------------------------------------
    // Function: batch_insert
    // Groups key-value pairs per rank and asynchronously inserts them.
    //------------------------------------------------------
    upcxx::future<> batch_insert(const std::vector<std::pair<Key, Value>> &entries) {
        std::unordered_map<int, std::vector<std::pair<Key, Value>>> grouped_inserts;
        for (const auto &entry : entries) {
            grouped_inserts[get_target_rank(entry.first)].push_back(entry);
            local_cache[entry.first] = entry.second;
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
        // Use when_all to aggregate all remote RPC futures,
        // then return a future that indicates completion.
        return upcxx::when_all(rpc_futures).then([](auto) { return upcxx::make_future(); });
    }

    //------------------------------------------------------
    // Function: insert
    // Convenience method to insert a single key-value pair.
    //------------------------------------------------------
    void insert(const Key &key, const Value &value) {
        batch_insert({{key, value}}).wait();
    }

    //------------------------------------------------------
    // Function: find
    // Retrieves a value given a key, using local cache when possible.
    //------------------------------------------------------
    bool find(const Key &key, Value &result) {
        auto cache_it = local_cache.find(key);
        if (cache_it != local_cache.end()) {
            result = cache_it->second;
            return true;
        }

        int target_rank = get_target_rank(key);
        if (target_rank == rank_id_) {
            auto it = local_map->find(key);
            if (it == local_map->end())
                return false;
            result = it->second;
            local_cache[key] = result;
            return true;
        } else {
            Value found_value = upcxx::rpc(target_rank,
                [](dist_map_t &lmap, const Key &key) -> Value {
                    auto it = lmap->find(key);
                    return (it != lmap->end()) ? it->second : Value{};
                }, local_map, key).wait();

            if (!(found_value == Value{})) {
                result = found_value;
                local_cache[key] = result;
                return true;
            }
            return false;
        }
    }

    //------------------------------------------------------
    // Function: process_requests
    // Advances progress on pending asynchronous UPC++ operations.
    //------------------------------------------------------
    void process_requests() {
        upcxx::progress(upcxx::progress_level::user);
    }
};
