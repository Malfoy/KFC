#include "strict_fstream.hpp"
#include "zstr.hpp"
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

#include "Kmers.cpp"
#include "Kmers.hpp"
#include "SuperKmerCount.cpp"
#include "SuperKmerCount.hpp"
#include "pow2.hpp"
#include "robin_hood.h"

using namespace std;

bool check = false;
robin_hood::unordered_map<string, uint64_t> real_count;
typedef robin_hood::unordered_flat_map<uint64_t,vector<SKC>> map;

uint64_t subminimizer_size(8);

uint64_t line_count(0);
bool aggressive_mode(false);

vector<omp_lock_t> MutexWall(4096);

Pow2<uint64_t> bucket_number(2 * subminimizer_size);

Pow2<uint64_t> offsetUpdateAnchor(2 * k);
Pow2<uint64_t> offsetUpdateAnchorMin(2 * minimizer_size);

//START LOW LEVEL FUNCTIONS

string getLineFasta(zstr::ifstream* in) {
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
	// cout<<(int)skc.weight<<"."<<(int)skc.size<<"	";
	string skmer(kmer2str(skc.sk, skc.size + 30));
	for (uint64_t i(0); i < skc.size; ++i) {
		// cout << getCanonical(skmer.substr(i, k)) << " " << (uint64_t)skc.counts[skc.size - 1 - i] << endl;
		if (check) {
			if (real_count[getCanonical(skmer.substr(i, k))] != (uint64_t)skc.counts[skc.size - 1 - i]) {
				cout << i << " " << skc.size - 1 - i << endl;
				cerr << "fail " << (skmer.substr(i, k)) << " " << revComp(skmer.substr(i, k)) << " real: " << real_count[getCanonical(skmer.substr(i, k))] << " output:"
				     << (uint64_t)skc.counts[skc.size - 1 - i] << endl;
				cin.get();
			}
		}
	}
}



void sort_buckets(vector<vector<SKC> >& buckets) {
	for (uint64_t i(0); i < buckets.size(); ++i) {
		sort(buckets[i].begin(), buckets[i].end(), [](const SKC& lhs, const SKC& rhs) {
			return lhs.weight > rhs.weight;
		});
	}
}



void dump_counting(map buckets[]) {
	// 	sort(buckets.begin(),buckets.end(),[ ]( const vector<SKC>& lhs, const vector<SKC>& rhs )
	// {
	//    return lhs.size()> rhs.size();
	// });
	for (uint64_t i(0); i < bucket_number.value(); ++i) {
		for (auto e:buckets[i]) {
			for (uint64_t ii(0); ii < e.second.size(); ++ii) {
				dump_count(e.second[ii]);
			}
		}
	}
	if (check) {
		cerr << "The results were checked" << endl;
	}
}



void dump_stats(map buckets[]) {
	uint64_t total_super_kmers(0);
	uint64_t total_kmers(0);
	uint64_t non_null_buckets(0);
	uint64_t null_buckets(0);
	uint64_t largest_bucket(0);

	//FOREACH BUCKETS
	for (uint64_t i(0); i < bucket_number.value(); ++i) {
		if (buckets[i].size() != 0) {
			for (auto e:buckets[i]) {
				largest_bucket = max(largest_bucket, (uint64_t)e.second.size());
				non_null_buckets++;
				total_super_kmers += e.second.size();
				for (uint64_t j(0); j < e.second.size(); ++j) {
					total_kmers += e.second[j].size;
				}
			}
		} else {
			null_buckets++;
		}
	}
	cout << endl;
	cout << "Empty buckets:	" << intToString(null_buckets) << endl;
	cout << "Useful buckets:	" << intToString(non_null_buckets) << endl;
	cout << "#Superkmer:	" << intToString(total_super_kmers) << endl;
	cout << "#Superkmer2:	" << intToString(nb_superkmer) << endl;
	cout << "#kmer:	" << intToString(total_kmers) << endl;
	cout << "super_kmer per useful buckets*1000:	" << intToString(total_super_kmers * 1000 / non_null_buckets) << endl;
	cout << "kmer per useful buckets*1000:	" << intToString(total_kmers * 1000 / non_null_buckets) << endl;
	cout << "kmer per super_kmer*1000:	" << intToString(total_kmers * 1000 / total_super_kmers) << endl;
	cout << "Largest_bucket:	" << intToString(largest_bucket) << endl;
}



void insert_kmers_into_bucket_last_chance(vector<kmer_full>& kmers, vector<SKC>& bucket, uint64_t minimizer, bool placed[], uint64_t bucket_offset) {
	// cout<<"go insert kmers"<<endl;
	uint64_t size_sk(kmers.size());
	uint64_t size_skc(bucket.size());

	//FOREACH KMER
	for (uint64_t ik = 0; ik < size_sk; ++ik) {
		if (not placed[ik]) {
			kmer_full& kmer = kmers[ik];
			// // FOREACH SUPERKMER
			for (uint64_t i = (bucket_offset); i < size_skc and (not placed[ik]); i++) {
				//Try to add the kmer into the superkmer counter
				placed[ik] = bucket[i].add_kmer(kmer);
			}
			// Create a new bucket if not placed
			if (not placed[ik]) {
				// FReeze the previous superker.
				if (size_skc > 0)
					bucket[size_skc - 1].close_compaction();
				// Read in the same strand than its canonical MINIMIZER
				if (minimizer == kmer.get_minimizer()) {
					bucket.push_back(SKC(kmer.kmer_s, kmer.get_minimizer_idx()));

				} else {
					// Reverse strand
					uint8_t start_idx = k - minimizer_size - kmer.get_minimizer_idx();
					bucket.push_back(SKC(kmer.kmer_rc, start_idx));
				}
				++size_skc;
			}
		}
	}

	bucket[size_skc - 1].close_compaction();
	kmers.clear();
}



void insert_kmers_into_bucket(vector<kmer_full>& kmers, vector<SKC>& bucket, uint64_t minimizer) {
	uint64_t size_sk(kmers.size());
	uint64_t size_skc(bucket.size());
	// vector<bool> placed(size_sk, false);
	bool placed[size_sk] = {false};
	uint64_t inserted    = 0;
	//FOREACH SUPERKMER
	for (uint64_t i = 0; i < size_skc; i++) {
		SKC& skc = bucket[i];
		// if(aggressive_mode and skc.weight<6){
		// 	break;
		// }
		// TODO: verify minimizer
		// uint64_t sk_mini = skc.get_minimizer();

		//FOREACH KMER
		for (uint64_t ik = 0; ik < size_sk; ++ik) {
			if (not placed[ik]) {
				kmer_full& kmer = kmers[ik];
				placed[ik]      = skc.add_kmer(kmer);
				if (placed[ik]) {
					inserted += 1;
				}
			}
		}

		if (i > 0) {
			if (bucket[i].weight >= bucket[i - 1].weight) {
				swap(bucket[i], bucket[i - 1]);
			}
		}
	}

	if (inserted == size_sk) {
		kmers.clear();
		return;
	}

	insert_kmers_into_bucket_last_chance(kmers, bucket, minimizer, placed, size_skc);
}



void count_line(const string& line, map buckets[]) {
	if (line.size() < k) {
		return;
	}
	vector<kmer_full> kmers;
	// Init Sequences
	uint64_t kmer_seq = (str2num(line.substr(0, k))), kmer_rc_seq(rcb(kmer_seq));
	uint64_t min_seq  = (str2num(line.substr(k - minimizer_size, minimizer_size))), min_rcseq(rcbc(min_seq, minimizer_size)), min_canon(min(min_seq, min_rcseq));
	// Init MINIMIZER
	int8_t relative_min_position;
	uint64_t minimizer = get_minimizer(kmer_seq, relative_min_position);
	bool multiple_min  = relative_min_position < 0;

	uint8_t position_minimizer_in_kmer;
	if (multiple_min)
		position_minimizer_in_kmer = (uint8_t)(-relative_min_position - 1);
	else
		position_minimizer_in_kmer = relative_min_position;

	uint64_t hash_mini = hash64shift(minimizer);
	kmers.push_back({relative_min_position, kmer_seq, kmer_rc_seq});
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
			insert_kmers_into_bucket(kmers, buckets[bucketindice][minimizer], minimizer);
			omp_unset_lock(&MutexWall[bucketindice % 4096]);

			minimizer                  = (min_canon);
			hash_mini                  = new_hash;
			position_minimizer_in_kmer = relative_min_position = 0;
			multiple_min                                       = false;
		}
		// duplicated MINIMIZER
		else if (new_hash == hash_mini) {
			multiple_min = true;
			position_minimizer_in_kmer += 1;
			relative_min_position = -((int8_t)position_minimizer_in_kmer) - 1;
		}
		//the previous MINIMIZER is outdated
		else if (position_minimizer_in_kmer >= k - minimizer_size) {
			uint64_t bucketindice = revhash(hash_mini) % bucket_number;
			omp_set_lock(&MutexWall[bucketindice % 4096]);
			insert_kmers_into_bucket(kmers, buckets[bucketindice][minimizer], minimizer);
			omp_unset_lock(&MutexWall[bucketindice % 4096]);
			// Search for the new MINIMIZER in the whole kmer
			minimizer    = get_minimizer(kmer_seq, relative_min_position);
			multiple_min = relative_min_position < 0;
			if (multiple_min)
				position_minimizer_in_kmer = (uint8_t)(-relative_min_position - 1);
			else
				position_minimizer_in_kmer = relative_min_position;

			hash_mini = hash64shift(minimizer);
		} else {
			position_minimizer_in_kmer++;
			if (multiple_min)
				relative_min_position--;
			else
				relative_min_position++;
		}

		// Normal add of the kmer into kmer list
		kmers.push_back({relative_min_position, kmer_seq, kmer_rc_seq});
	}

	uint64_t bucketindice = revhash(hash_mini) % bucket_number;
	omp_set_lock(&MutexWall[bucketindice % 4096]);
	insert_kmers_into_bucket(kmers, buckets[bucketindice][minimizer], minimizer);
	omp_unset_lock(&MutexWall[bucketindice % 4096]);
}



void read_fasta_file(const string& filename, map buckets[]) {
	for (uint64_t i(0); i < 4096; ++i) {
		omp_init_lock(&MutexWall[i]);
	}
	zstr::ifstream in(filename);
	uint8_t nb_core(8);
	if (check) {
		nb_core = 1;
	}
#pragma omp parallel num_threads(nb_core)
	{
		string line;
		while (in.good()) {
#pragma omp critical(input)
			{
				line = getLineFasta(&in);
			}
			count_line(line, buckets);
			line_count++;
			// #pragma omp critical (aggro)
			// {
			// 	if(not aggressive_mode and nb_superkmer>(1000*10*1)){
			// 		cout<<"AGGRO MODE"<<endl;
			// 		aggressive_mode=true;
			// 		sort_buckets(buckets);
			// 	}
			// }

			if (line_count % 100000 == 0) {
				// cout<<nb_kmer_read<<" "<<nb_superkmer<<" "<<(double)nb_kmer_read/nb_superkmer<<endl,
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
	// if (argc > 3) {
	// 	minimizer_size = stoi(argv[3]);
	// }

	if (mode > 1) {
		check = true;
		cout << "LETS CHECK THE RESULTS" << endl;
	}

	map buckets[bucket_number.value()];
	cout << "\n\n\nI count " << argv[1] << endl;
	cout << "Minimizer size:	" << minimizer_size << endl;
	cout << "Number of bucket:	" << bucket_number.value() << endl;

	read_fasta_file(argv[1], buckets);
	cout << endl;
	if (mode % 2 == 0) {
		dump_counting(buckets);
	}
	dump_stats(buckets);
	exit(0);
}
