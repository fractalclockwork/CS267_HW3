#pragma once
#include <upcxx/upcxx.hpp>
#include <unordered_map>
#include <string>
#include "kmer_t.hpp"

class DistributedHashMap {
private:
    using dobj_map_t = upcxx::dist_object<std::unordered_map<std::string, kmer_pair>>;
    dobj_map_t local_map;
    size_t table_size_;
    int rank_id_;
    int world_size_;

    int get_target_rank(const std::string &key) const {
        return std::hash<std::string>{}(key) % world_size_;
    }

public:
    DistributedHashMap(size_t table_size, int rank_id, int world_size)
        : table_size_(table_size), rank_id_(rank_id), world_size_(world_size), local_map({}) {}

    upcxx::future<> insert(const std::string &key, const kmer_pair &value) {
        int target_rank = get_target_rank(key);
        if (target_rank == rank_id_) {
            local_map->insert({key, value});
            return upcxx::make_future();
        } else {
            return upcxx::rpc(target_rank,
                [](dobj_map_t &lmap, const std::string &key, const kmer_pair &value) {
                    lmap->insert({key, value});
                },
                local_map, key, value);
        }
    }

    upcxx::future<bool> find(const std::string &key, kmer_pair &result) {
        int target_rank = get_target_rank(key);
        if (target_rank == rank_id_) {
            auto it = local_map->find(key);
            if (it == local_map->end()) {
                return upcxx::make_future(false);
            }
            result = it->second;
            return upcxx::make_future(true);
        } else {
            return upcxx::rpc(target_rank,
                [](dobj_map_t &lmap, const std::string &key) -> kmer_pair {
                    auto it = lmap->find(key);
                    return (it != lmap->end()) ? it->second : kmer_pair();
                },
                local_map, key
            ).then([&result](kmer_pair found_kmer) {
                if (!found_kmer.kmer_str().empty()) {
                    result = found_kmer;
                    return true;
                }
                return false;
            });
        }
    }

    void process_requests() {
        upcxx::progress(upcxx::progress_level::user);
    }

    void synchronize() {
        upcxx::barrier();
    }
};