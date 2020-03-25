#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <fstream>
#include <gatbl/common.hpp>
#include <gatbl/kmer.hpp>
#include <iostream>
#include <math.h>
#include <mutex>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <tmmintrin.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>

#include "Kmers.hpp"
#include "SuperKmerCount.hpp"
#include "pow2.hpp"
#include "robin_hood.h"

using namespace std;

bool check = false;
robin_hood::unordered_map<string, uint64_t> real_count;

uint64_t subminimizer_size(minimizer_size - 1);
uint64_t nb_kmer_read(0);
uint64_t line_count(0);

vector<omp_lock_t> MutexWall(4096);

Pow2<uint64_t> bucket_number(2 * subminimizer_size);

Pow2<uint64_t> offsetUpdateAnchor(2 * k);
Pow2<uint64_t> offsetUpdateAnchorMin(2 * minimizer_size);

//START LOW LEVEL FUNCTIONS

string getLineFasta(ifstream* in) {
	string line, result;
	getline(*in, line);
	char c = static_cast<char>(in->peek());
	while (c != '>' and c != EOF) {
		getline(*in, line);
		result += line;
		c = static_cast<char>(in->peek());
	}
	return result;
}

uint32_t revhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x2c1b3c6d;
	x = ((x >> 16) ^ x) * 0x297a2d39;
	x = ((x >> 16) ^ x);
	return x;
}

uint32_t unrevhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x0cf0b109; // PowerMod[0x297a2d39, -1, 2^32]
	x = ((x >> 16) ^ x) * 0x64ea2d65;
	x = ((x >> 16) ^ x);
	return x;
}

uint64_t revhash(uint64_t x) {
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x);
	return x;
}

uint64_t unrevhash(uint64_t x) {
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
}

inline uint64_t nuc2int(char c) {
	return (c / 2) % 4;
}

inline uint64_t nuc2intrc(char c) {
	return ((c / 2) % 4) ^ 2;
}

inline void updateK(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchor;
}

inline void add_nuc_superkmer(SKC& min, char nuc) {
	min.sk <<= 2;
	min.sk += nuc2int(nuc);
	min.counts[min.size++] = 0;
}

inline void updateM(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchorMin;
}

inline void updateRCK(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * k - 2));
}

inline void updateRCM(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * minimizer_size - 2));
}

char revCompChar(char c) {
	switch (c) {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
	}
	return 'A';
}

string revComp(const string& s) {
	string rc(s.size(), 0);
	for (int i((int)s.length() - 1); i >= 0; i--) {
		rc[s.size() - 1 - i] = revCompChar(s[i]);
	}
	return rc;
}

string getCanonical(const string& str) {
	return (min(str, revComp(str)));
}

string intToString(uint64_t n) {
	if (n < 1000) {
		return to_string(n);
	}
	string end(to_string(n % 1000));
	if (end.size() == 3) {
		return intToString(n / 1000) + "," + end;
	}
	if (end.size() == 2) {
		return intToString(n / 1000) + ",0" + end;
	}
	return intToString(n / 1000) + ",00" + end;
}

//START HIGH LEVEL FUNCTIONS

void dump_count(const SKC& skc) {
	cout<<(int)skc.weight<<"."<<(int)skc.size<<"	";
	string skmer(kmer2str(skc.sk, skc.size + 30));
	for (uint64_t i(0); i < skc.size; ++i) {
		// cout << getCanonical(skmer.substr(i, k)) << " " << (uint64_t)skc.counts[skc.size - 1 - i] << endl;
		if (check) {
			// if (real_count[getCanonical(skmer.substr(i, k))] != (uint64_t)skc.counts[skc.size - 1 - i]) {
			// 	cout << i << " " << skc.size - 1 - i << endl;
			// 	cerr << "fail " << (skmer.substr(i, k)) << " " << revComp(skmer.substr(i, k)) << " real: " << real_count[getCanonical(skmer.substr(i, k))] << " output:"
			// 	     << (uint64_t)skc.counts[skc.size - 1 - i] << endl;
			// 	cin.get();
			// }
		}
	}
}

void dump_counting(const vector<vector<SKC> >& buckets) {
	auto copy(buckets);
	sort(copy.begin(),copy.end(),[ ]( const vector<SKC>& lhs, const vector<SKC>& rhs )
{
   return lhs.size()> rhs.size();
});
	for (uint64_t i(0); i < copy.size(); ++i) {
		if(copy[i].size()!=0){
			cout<<endl<<"go buckets"<<endl;
		}
		for (uint64_t ii(0); ii < copy[i].size(); ++ii) {
			dump_count(copy[i][ii]);
		}
		if(copy[i].size()!=0){
			cin.get();
		}
	}
	if (check) {
		cerr << "The results were checked" << endl;
	}
}

void dump_stats(const vector<vector<SKC> >& buckets) {
	uint64_t total_super_kmers(0);
	uint64_t total_kmers(0);
	uint64_t non_null_buckets(0);
	uint64_t null_buckets(0);
	uint64_t largest_bucket(0);

	//FOREACH BUCKETS
	for (uint64_t i(0); i < buckets.size(); ++i) {
		if (buckets[i].size() != 0) {
			largest_bucket = max(largest_bucket, (uint64_t)buckets[i].size());
			non_null_buckets++;
			total_super_kmers += buckets[i].size();
			for (uint64_t j(0); j < buckets[i].size(); ++j) {
				total_kmers += buckets[i][j].size;
			}
		} else {
			null_buckets++;
		}
	}
	cout << endl;
	cout << "Empty buckets:	" << intToString(null_buckets) << endl;
	cout << "Useful buckets:	" << intToString(non_null_buckets) << endl;
	cout << "#Superkmer:	" << intToString(total_super_kmers) << endl;
	cout << "#kmer:	" << intToString(total_kmers) << endl;
	cout << "super_kmer per useful buckets:	" << intToString(total_super_kmers / non_null_buckets) << endl;
	cout << "kmer per useful buckets:	" << intToString(total_kmers / non_null_buckets) << endl;
	cout << "kmer per super_kmer:	" << intToString(total_kmers / total_super_kmers) << endl;
	cout << "Largest_bucket:	" << intToString(largest_bucket) << endl;
}


int compact(SKC& super_kmer, kmer_full kmer) {
	cout << "COMPACT" << endl;
	if (super_kmer.size >= 13) {
		return -1;
	}
	uint64_t superkmer = super_kmer.sk;
	uint64_t end_sk    = superkmer;
	end_sk &= (((uint64_t)1 << 60) - 1);
	uint64_t beg_k = kmer.kmer_s >> 2;

	print_kmer(end_sk, 30);
	print_kmer(kmer.kmer_s, 31);
	print_kmer(kmer.kmer_rc, 31);
	if (end_sk == beg_k) {
		super_kmer.sk <<= 2;
		super_kmer.sk += (((kmer.kmer_s % 4)));
		super_kmer.counts[super_kmer.size++] = 1;
		super_kmer.minimizer_idx++;
		return super_kmer.size;
	}
	beg_k = kmer.kmer_rc >> 2;
	if (end_sk == beg_k) {
		super_kmer.sk <<= 2; //TODO IDONOTKNOW
		super_kmer.sk += (((kmer.kmer_rc % 4)));
		super_kmer.counts[super_kmer.size++] = 1;
		return super_kmer.size;
	}
	cout << "COMPACT FAIL !!!!" << endl;
	cin.get();
	return -1;
}

void insert_kmers_into_bucket_last_chance(vector<kmer_full>& kmers, vector<SKC>& bucket, uint64_t minimizer, vector<bool> placed = vector<bool>(0, false)) {
	// cout<<"go insert kmers"<<endl;
	uint64_t size_sk(kmers.size());
	uint64_t size_skc(bucket.size());
	if (placed.size() != size_sk)
		placed = vector<bool>(size_sk, false);

	//FOREACH KMER
	for (uint64_t ik = 0; ik < size_sk; ++ik) {
		kmer_full& kmer = kmers[ik];
		// // FOREACH SUPERKMER
		for (uint64_t i = (0); i < size_skc and (not placed[ik]); i++) {
			//Try to add the kmer into the superkmer counter
			placed[ik] = bucket[i].add_kmer(kmer);
		}
		// Create a new bucket if not placed
		if (not placed[ik]) {
			// Read in the same strand than its canonical MINIMIZER
			if (minimizer == kmer.get_minimizer()) {
				bucket.push_back(SKC(kmer.kmer_s, kmer.minimizer_idx));
			}else {
				// Reverse strand
				uint8_t start_idx = k - minimizer_size - kmer.minimizer_idx;
				bucket.push_back(SKC(kmer.kmer_rc, start_idx));
			}
			size_skc++;
		}
	}
	kmers.clear();
}


void insert_kmers_into_bucket(vector<kmer_full>& kmers, vector<SKC>& bucket, uint64_t minimizer) {
	uint64_t size_sk(kmers.size());
	uint64_t size_skc(bucket.size());
	vector<bool> placed(size_sk, false);
	uint64_t inserted = 0;

	//FOREACH SUPERKMER
	for (uint64_t i = 0; i < size_skc; i++) {
		SKC& skc = bucket[i];

		//FOREACH KMER
		for (uint64_t ik = 0; ik < size_sk; ++ik) {
			kmer_full& kmer = kmers[ik];

			if (not placed[ik]) {
				placed[ik] = skc.add_kmer(kmer);
				if (placed[ik]){
					inserted += 1;

				}

			}
		}
		if(i>0){
			if(bucket[i].weight>=1.5*bucket[i-1].weight){
				swap(bucket[i],bucket[i-1]);
			}
		}
	}

	if (inserted == size_sk) {
		kmers.clear();
		return;
	}

	vector<SKC> new_skc = vector<SKC>();
	insert_kmers_into_bucket_last_chance(kmers, new_skc, minimizer, placed);
	bucket.insert(bucket.end(), new_skc.begin(), new_skc.end());
}


void count_line(const string& line, vector<vector<SKC> >& buckets) {
	if (line.size() < k) {
		return;
	}
	vector<kmer_full> kmers;
	// Init Sequences
	uint64_t kmer_seq = (str2num(line.substr(0, k))), kmer_rc_seq(rcb(kmer_seq));
	uint64_t min_seq  = (str2num(line.substr(k - minimizer_size, minimizer_size))), min_rcseq(rcbc(min_seq, minimizer_size)), min_canon(min(min_seq, min_rcseq));
	// Init MINIMIZER
	uint64_t position_min;
	uint64_t position_minimizer_in_kmer;
	uint64_t minimizer         = get_minimizer(kmer_seq, position_min);
	position_minimizer_in_kmer = position_min;

	uint64_t hash_mini = hash64shift(minimizer);
	kmers.push_back({static_cast<uint8_t>(position_minimizer_in_kmer), kmer_seq, kmer_rc_seq});
	if (check) {
		real_count[getCanonical(line.substr(0, k))]++;
	}
	uint64_t line_size = line.size();
	for (uint64_t i = 0; i + k < line_size; ++i) {
		if (check) {
			real_count[getCanonical(line.substr(i + 1, k))]++;
		}
		// Update KMER and MINIMIZER candidate with the new letter
		updateK(kmer_seq, line[i + k]);
		updateRCK(kmer_rc_seq, line[i + k]);
		updateM(min_seq, line[i + k]);
		updateRCM(min_rcseq, line[i + k]);
		min_canon = (min(min_seq, min_rcseq));

		//THE NEW mmer is a MINIMIZER
		uint64_t new_hash = (hash64shift(min_canon));
		if (new_hash < hash_mini) {
			uint64_t bucketindice = revhash(hash_mini) % bucket_number;
			omp_set_lock(&MutexWall[bucketindice % 4096]);
			insert_kmers_into_bucket(kmers, buckets[bucketindice], minimizer);
			omp_unset_lock(&MutexWall[bucketindice % 4096]);

			minimizer                  = (min_canon);
			hash_mini                  = new_hash;
			position_min               = i + k - minimizer_size + 1;
			position_minimizer_in_kmer = 0;
		} else {
			//the previous MINIMIZER is outdated
			if (position_minimizer_in_kmer>=k-minimizer_size) {
				uint64_t bucketindice = revhash(hash_mini) % bucket_number;
				omp_set_lock(&MutexWall[bucketindice % 4096]);
				insert_kmers_into_bucket(kmers, buckets[bucketindice], minimizer);
				omp_unset_lock(&MutexWall[bucketindice % 4096]);
				// Search for the new MINIMIZER in the whole kmer
				minimizer = get_minimizer(kmer_seq, position_min);

				position_minimizer_in_kmer = position_min;
				hash_mini                  = hash64shift(minimizer);
				position_min += i + 1;
			}else{
				position_minimizer_in_kmer++;
			}
		}
		// Normal add of the kmer into kmer list
		kmers.push_back({static_cast<uint8_t>(position_minimizer_in_kmer), kmer_seq, kmer_rc_seq});
	}

	uint64_t bucketindice = revhash(hash_mini) % bucket_number;
	omp_set_lock(&MutexWall[bucketindice % 4096]);
	insert_kmers_into_bucket(kmers, buckets[bucketindice], minimizer);
	omp_unset_lock(&MutexWall[bucketindice % 4096]);
}

void read_fasta_file(const string& filename, vector<vector<SKC> >& buckets) {
	for (uint64_t i(0); i < 4096; ++i) {
		omp_init_lock(&MutexWall[i]);
	}
	ifstream in(filename);
	uint8_t nb_core(4);
	if (check) {
		nb_core = 1;
	}
#pragma omp parallel num_threads(nb_core)
	{
		string line;
		while (in.good()) {
#pragma omp critical
			{
				line = getLineFasta(&in);
			}
			count_line(line, buckets);
			line_count++;
			if (line_count % 100000 == 0) {
				cerr << "-" << flush;
			}
		}
	}
}

int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "[fasta file]" << endl;
		cout << "0 DEFAULT, OUTPUT THE COUNT WITH NO CHECKING" << endl;
		cout << "1 PERFORMANCE MODE NO COUNT OUTPUT" << endl;
		cout << "2 DEBUG MODE COUNT AND CHECK WITH HASHTABLE" << endl;
		exit(0);
	}
	int mode(0);
	if (argc > 2) {
		mode = stoi(argv[2]);
	}

	if (mode > 1) {
		check = true;
		cout<<"LETS CHECK THE RESULTS"<<endl;
	}
	vector<vector<SKC> > buckets(bucket_number.value());
	read_fasta_file(argv[1], buckets);
	cout << endl;
	if (mode % 2 == 0) {
		dump_counting(buckets);
	}
	dump_stats(buckets);
	exit(0);
}
