#include <upcxx/upcxx.hpp>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>

//======================================================
// DistributedHash: Low-level distributed hash table
//======================================================
template<typename Key, typename Value>
class DistributedHash {
public:
    using HashTable = std::unordered_map<Key, Value>;
    using DistMap = upcxx::dist_object<HashTable>;

private:
    DistMap local_data;
    int rank_id_;
    int world_size_;

public:
    DistributedHash(int rank_id, int world_size)
        : rank_id_(rank_id), world_size_(world_size), local_data(HashTable{}) {}

    // Insert key/value pair (local or remote)
    void insert(const Key& key, const Value& value) {
        int target = get_target_rank(key);
        if (target == rank_id_) {
            local_data->emplace(key, value);
        } else {
            upcxx::rpc_ff(target,
                [](upcxx::dist_object<HashTable>& dobj, const Key& key, const Value& value) {
                    dobj->emplace(key, value);
                },
                local_data, key, value);
        }
    }

    // Lookup operation
    upcxx::future<Value> find(const Key& key) {
        int target = get_target_rank(key);
        if (target == rank_id_) {
            auto it = local_data->find(key);
            return upcxx::make_future(it != local_data->end() ? it->second : Value{});
        } else {
            return upcxx::rpc(target,
                [](upcxx::dist_object<HashTable>& dobj, const Key& key) -> Value {
                    auto it = dobj->find(key);
                    return (it != dobj->end()) ? it->second : Value{};
                },
                local_data, key);
        }
    }

    // Process pending UPC++ RPCs
    void process_requests() {
        upcxx::progress(upcxx::progress_level::user);
    }

private:
    int get_target_rank(const Key& key) const {
        return std::hash<Key>{}(key) % world_size_;
    }
};

//======================================================
// DistributedHashMap: Higher-level abstraction
//======================================================
template<typename Key, typename Value>
class DistributedHashMap {
private:
    upcxx::dist_object<DistributedHash<Key, Value>> d_hash;
    std::unordered_map<Key, Value> local_cache;
    int rank_id_;
    int world_size_;

    int get_target_rank(const Key& key) const {
        return std::hash<Key>{}(key) % world_size_;
    }

public:
    DistributedHashMap(int rank_id, int world_size)
        : rank_id_(rank_id), world_size_(world_size), d_hash(DistributedHash<Key, Value>(rank_id, world_size)) {}

    // Batch Insert: Efficient grouped insertions using UPC++ futures
    upcxx::future<> batch_insert(const std::vector<std::pair<Key, Value>>& entries) {
        std::unordered_map<int, std::vector<std::pair<Key, Value>>> grouped_batches;
        for (const auto& entry : entries) {
            grouped_batches[get_target_rank(entry.first)].push_back(entry);
            local_cache[entry.first] = entry.second; // Local caching
        }

        std::vector<upcxx::future<>> rpc_futures;

        for (auto& [target, batch] : grouped_batches) {
            if (target == rank_id_) {
                for (const auto& kv : batch) {
                    d_hash->insert(kv.first, kv.second);
                }
            } else {
                rpc_futures.push_back(
                    upcxx::rpc(target,
                        [](upcxx::dist_object<DistributedHash<Key, Value>>& remote_hash,
                           std::vector<std::pair<Key, Value>> batch) {
                            for (const auto& kv : batch) {
                                remote_hash->insert(kv.first, kv.second);
                            }
                        },
                        d_hash, std::move(batch)));
            }
        }

        if (rpc_futures.empty()) {
            return upcxx::make_future();
        } else {
            return upcxx::when_all(std::move(rpc_futures)).then([](auto) { return upcxx::make_future(); });
        }
    }

    // Single Insert Wrapper
    void insert(const Key& key, const Value& value) {
        batch_insert({{key, value}}).wait();
    }

    // Lookup Operation
    upcxx::future<Value> find(const Key& key) {
        auto cache_it = local_cache.find(key);
        if (cache_it != local_cache.end()) {
            return upcxx::make_future(cache_it->second);
        }
        return d_hash->find(key).then([this, key](Value result) {
            if (!(result == Value{})) {
                local_cache[key] = result;
            }
            return result;
        });
    }

    // Process pending UPC++ requests
    void process_requests() {
        d_hash->process_requests();
    }
};

//======================================================
// Main: Testing the Optimized DHM
//======================================================
int main() {
    upcxx::init();

    int rank = upcxx::rank_me();
    int world_size = upcxx::rank_n();

    DistributedHashMap<int, std::string> dhm(rank, world_size);

    // Generate test data (100 key-value pairs)
    std::vector<std::pair<int, std::string>> entries;
    for (int i = 0; i < 100; i++) {
        entries.push_back({i, "Value_" + std::to_string(i)});
    }

    // Perform batch insertion
    dhm.batch_insert(entries).wait();

    upcxx::barrier();

    // Lookup and print results on rank 0
    if (rank == 0) {
        for (int i = 0; i < 100; i++) {
            std::string value = dhm.find(i).wait();
            std::cout << "Key: " << i << " Value: " << value << std::endl;
        }
    }

    upcxx::finalize();
    return 0;
}
