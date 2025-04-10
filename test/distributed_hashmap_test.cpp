#include <upcxx/upcxx.hpp>
#include <iostream>
#include <unordered_map>
#include <string>

//--------------------------------------------------------
// Generic Distributed Hash Map using UPC++ PGAS Model
//--------------------------------------------------------
template<typename Key, typename Value, typename Hash = std::hash<Key>>
class DistributedHashMap {
private:
    // Local hash table stored on each rank.
    using LocalMap = std::unordered_map<Key, Value, Hash>;
    upcxx::dist_object<LocalMap> local_data;

    // Store rank and world size for steering.
    const int rank;
    const int world_size;

    // Compute the target rank for a key.
    int get_target_rank(const Key& key) const {
        return Hash{}(key) % world_size;
    }

public:
    // Default constructor initializes rank and world size automatically.
    DistributedHashMap()
        : rank(upcxx::rank_me()),
          world_size(upcxx::rank_n()),
          local_data(LocalMap{}) {}

    // Asynchronous put operation (fire-and-forget)
    void put(const Key& key, const Value& value) {
        int target = get_target_rank(key);
        if (target == rank) {
            local_data->insert_or_assign(key, value);
        } else {
            upcxx::rpc_ff(target,
                // Remote lambda: insert or update the key/value pair.
                [](upcxx::dist_object<LocalMap>& remote_data,
                   const Key& key, const Value& value) {
                    remote_data->insert_or_assign(key, value);
                },
                local_data, key, value);
        }
    }

    // Asynchronous get operation: returns a future holding the Value.
    upcxx::future<Value> get(const Key& key) const {
        int target = get_target_rank(key);
        if (target == rank) {
            auto it = local_data->find(key);
            if (it != local_data->end())
                return upcxx::make_future(it->second);
            else
                return upcxx::make_future(Value{});  // Return default-constructed Value
        } else {
            return upcxx::rpc(target,
                // Remote lambda: lookup key and return result.
                [](upcxx::dist_object<LocalMap>& remote_data,
                   const Key& key) -> Value {
                    auto it = remote_data->find(key);
                    return (it != remote_data->end()) ? it->second : Value{};
                },
                local_data, key);
        }
    }

    // Optionally expose progress to help drive UPC++ communication.
    inline void progress() { 
        upcxx::progress(upcxx::progress_level::user); 
    }
};


//--------------------------------------------------------
// Testing the Generic Distributed Hash Map
//--------------------------------------------------------
int main() {
    upcxx::init();

    {
        // Create a distributed hash map with int keys and std::string values.
        DistributedHashMap<int, std::string> dhm;

        // Insert 100 key/value pairs.
        for (int i = 0; i < 100; ++i) {
            dhm.put(i, "Value_" + std::to_string(i));
        }

        // Synchronize to ensure all RPCs have been issued.
        upcxx::barrier();

        // On rank 0, lookup each key and print the value.
        if (upcxx::rank_me() == 0) {
            for (int i = 0; i < 100; ++i) {
                std::string value = dhm.get(i).wait();
                std::cout << "Key: " << i << " Value: " << value << std::endl;
            }
        }

        upcxx::barrier();
    }

    upcxx::finalize();
    return 0;
}
