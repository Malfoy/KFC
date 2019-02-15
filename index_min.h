#ifndef __INDEX8MIN
#define __INDEX8MIN



#include <stdio.h>
#include <fstream>
#include <vector>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <algorithm>
#include "BBHash/BooPHF.h"


using namespace std;



typedef boomphf::SingleHashFunctor<uint64_t>  hasher;
typedef boomphf::mphf<  uint64_t, hasher > MPHF;






struct value{
	uint64_t kmer;
	uint8_t count;
};





class index_full{
public:
	uint32_t kmer_size;
	MPHF Hash;
	vector<value> Values;

	index_full(const vector<uint64_t>& V){
		Values.resize(V.size());
		auto data_iterator = boomphf::range(static_cast<const uint64_t*>(&((V)[0])), static_cast<const uint64_t*>((&(V)[0])+V.size()));
		Hash= boomphf::mphf<uint64_t,hasher>(V.size(),data_iterator,4,5,false);
	}
	void insert(uint64_t);
	void dump_counting();
	void insert_seq(const string&  read);

};




class index_min{
public:
	uint kmer_size;
	uint minimizer_size;

	vector<index_full> Index;

	index_min(vector<uint64_t>& V);
	void insert(uint64_t kmer);
	void dump_counting();

};



#endif
