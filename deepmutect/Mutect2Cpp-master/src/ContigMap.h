//
// Created by hlf on 6/29/22.
//

#ifndef MUTECT2CPP_MASTER_CONTIGMAP_H
#define MUTECT2CPP_MASTER_CONTIGMAP_H

#include "parallel_hashmap/phmap.h"
#include <vector>
#include <string>
#include <stdexcept>
#include <xxhash.hpp>

struct hash_contig {
	xxh::hash64_t operator()(const std::string &s1) const;
};

class ContigMap {
private:
	/**
	 * intToString: {"", "chr1", "chr2", ...}
	 * stringToInt: {"", -1}, {"chr1", 0}, {"chr2", 1}, ...
	 */
	static std::vector<std::string> intToString;
	static phmap::flat_hash_map<std::string, int, hash_contig> stringToInt;
	static int mapSize;

public:
	static void initial(int reserveSize);

	static int getContigInt(const std::string &key);

	static std::string getContigString(int key);

	static void insertPair(int intKey, const std::string &stringKey);
};

#endif //MUTECT2CPP_MASTER_CONTIGMAP_H
