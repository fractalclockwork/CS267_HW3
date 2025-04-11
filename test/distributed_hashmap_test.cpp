#include <upcxx/upcxx.hpp>
#include <iostream>
#include <vector>
#include <string>
#include <functional>
#include <stdexcept>

//--------------------------------------------------------------------
// VectorHashTable:
// A simple open-addressing hash table implemented over std::vector.
//--------------------------------------------------------------------
template<typename Key, typename Value, typename Hash = std::hash<Key>>
class VectorHashTable {
public:
    // Each table slot is an Entry that tracks occupancy.
    struct Entry {
        bool occupied;
        Key key;
        Value value;
        Entry() : occupied(false), key(), value() {}
    };

private:
    std::vector<Entry> table;
    size_t capacity;
    Hash hash_fn;

public:
    // Construct with a fixed capacity (default is 256).
    explicit VectorHashTable(size_t cap = 256)
        : capacity(cap), table(cap) {}

    // Insert or assign a key/value pair.
    void insert_or_assign(const Key& key, const Value& value) {
        size_t start = hash_fn(key) % capacity;
        for (size_t i = 0; i < capacity; i++) {
            size_t idx = (start + i) % capacity;
            if (!table[idx].occupied) { 
                table[idx].occupied = true;
                table[idx].key = key;
                table[idx].value = value;
                return;
            } else if (table[idx].occupied && table[idx].key == key) {
                table[idx].value = value;
                return;
            }
        }
        // If we loop through the entire table, it is full.
        throw std::overflow_error("Hash table is full");
    }

    // Find a key; returns pointer to value if found, or nullptr otherwise.
    Value* find(const Key& key) {
        size_t start = hash_fn(key) % capacity;
        for (size_t i = 0; i < capacity; i++) {
            size_t idx = (start + i) % capacity;
            if (table[idx].occupied && table[idx].key == key) {
                return &table[idx].value;
            } else if (!table[idx].occupied) {
                // Slot is empty; key does not exist.
                return nullptr;
            }
        }
        return nullptr;
    }

    // Const version of find.
    const Value* find(const Key& key) const {
        size_t start = hash_fn(key) % capacity;
        for (size_t i = 0; i < capacity; i++) {
            size_t idx = (start + i) % capacity;
            if (table[idx].occupied && table[idx].key == key) {
                return &table[idx].value;
            } else if (!table[idx].occupied) {
                return nullptr;
            }
        }
        return nullptr;
    }
};

//--------------------------------------------------------------------
// DistributedHashMap:
// A generic distributed hash map that uses VectorHashTable as local
// storage. The UPC++ PGAS model is used to distribute the data across
// ranks via upcxx::dist_object.
//--------------------------------------------------------------------
template<typename Key, typename Value, typename Hash = std::hash<Key>>
class DistributedHashMap {
private:
    using LocalTable = VectorHashTable<Key, Value, Hash>;
    upcxx::dist_object<LocalTable> local_data;
    const int rank;
    const int world_size;

    // Compute the rank that should own the key.
    int get_target_rank(const Key& key) const {
        return Hash{}(key) % world_size;
    }

public:
    // Constructor initializes the local table with a fixed capacity.
    DistributedHashMap()
        : rank(upcxx::rank_me()),
          world_size(upcxx::rank_n()),
          local_data(LocalTable(256)) {}

    // Asynchronous put: inserts or updates a key/value pair.
    void put(const Key& key, const Value& value) {
        int target = get_target_rank(key);
        if (target == rank) {
            local_data->insert_or_assign(key, value);
        } else {
            upcxx::rpc_ff(target,
                // Remote lambda to perform the insertion.
                [](upcxx::dist_object<LocalTable>& remote_table,
                   const Key& key, const Value& value) {
                    remote_table->insert_or_assign(key, value);
                },
                local_data, key, value);
        }
    }

    // Asynchronous get: returns a UPC++ future to the value.
    upcxx::future<Value> get(const Key& key) const {
        int target = get_target_rank(key);
        if (target == rank) {
            Value* result = local_data->find(key);
            if (result) {
                return upcxx::make_future(*result);
            } else {
                return upcxx::make_future(Value{}); // Default value if not found.
            }
        } else {
            return upcxx::rpc(target,
                // Remote lambda to perform the lookup.
                [](upcxx::dist_object<LocalTable>& remote_table,
                   const Key& key) -> Value {
                       const Value* result = remote_table->find(key);
                       return result ? *result : Value{};
                },
                local_data, key);
        }
    }

    // Optional progress method to advance UPC++ communication.
    inline void progress() {
        upcxx::progress(upcxx::progress_level::user);
    }
};

//--------------------------------------------------------------------
// main(): Testing the Generic Distributed Hash Map with Vector
//--------------------------------------------------------------------
int main() {
    upcxx::init();

    {
        // Create a distributed hash map using int keys and std::string values.
        DistributedHashMap<int, std::string> dhm;

        // Each rank inserts 100 key/value pairs.
        for (int i = 0; i < 100; ++i) {
            dhm.put(i, "Value_" + std::to_string(i));
        }

        // Synchronize to ensure all remote invocations complete.
        upcxx::barrier();

        // On rank 0, perform lookups for keys 0 through 99 and print results.
        if (upcxx::rank_me() == 0) {
            for (int i = 0; i < 100; ++i) {
                std::string result = dhm.get(i).wait();
                std::cout << "Key: " << i << " Value: " << result << std::endl;
            }
        }

        upcxx::barrier();
    }

    upcxx::finalize();
    return 0;
}
