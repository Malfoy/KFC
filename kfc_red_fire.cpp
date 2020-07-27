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
#include "SuperKmerCount.cpp"
#include "buckets.cpp"
#include "pow2.hpp"
#include "robin_hood.h"
#include "sparse_map.h"



using namespace std;



const uint64_t subminimizer_size(15);//min 5 with 21
uint64_t line_count(0);
bool aggressive_mode(false);
vector<omp_lock_t> MutexWall(4096);
const Pow2<uint64_t> bucket_number(2 * subminimizer_size);
Pow2<uint64_t> offsetUpdateAnchor(2 * k);
const Pow2<uint64_t> offsetUpdateAnchorMin(2 * minimizer_size);
uint16_t abundance_mini[1<<(2*minimizer_size)];
vector<Bucket> bucket_menus[1<<(2*subminimizer_size)];
uint64_t nb_core(8);



//START LOW LEVEL FUNCTIONS
uint64_t hash64shift(uint64_t key) {
	// uint64_t ab=abundance_mini[key];
	uint64_t ab=0;
	// ab<<=32;

	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return ab+(uint32_t)key;
}


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



void dump_count_bucket(const vector<SKC>& v,uint64_t mini){
	string out,skm,prefix,suffix,minimizer(kmer2str(mini,minimizer_size));
	// FOREACH SUPERKMER
	for(uint isk(0);isk<v.size();++isk){
		SKC skc(v[isk]);
		skm=kmer2str(skc.sk,30+skc.size-minimizer_size);
		prefix=skm.substr(0,skm.size()-skc.minimizer_idx);
		suffix=skm.substr(prefix.size());
		skm=prefix+minimizer+suffix;
		//FOREACH KMER WITHIN THE SUPERKMER
		for (uint64_t i(0); i < skc.size; ++i) {
			out+=skm.substr(i,31)+' '+to_string(skc.counts[i])+'\n';
			if(check){
				if(real_count[getCanonical(skm.substr(i,31))]!=skc.counts[i]){
					cout<<skm.substr(i,31)<<" "<<to_string(skc.counts[i]);
					cout<<"	instead of ";
					cout<<real_count[getCanonical(skm.substr(i,31))]<<endl;
					counting_errors++;
				}else{
					real_count[getCanonical(skm.substr(i,31))]=0;
				}
			}
		}
	}
	#pragma omp critical (outcerr)
	{
		// cerr<<out;
	}
}



void dump_counting() {
	// #pragma omp parallel for
	for (uint64_t i=0; i < bucket_number.value(); ++i) {
		string toprint;
		for(auto e:bucket_menus[i]){
			e.print_kmers(toprint,kmer2str(i,minimizer_size));
		}
	}
	if(check){
		for (auto e:real_count) {
			if(e.second!=0){
				cout<<"I forgot	"<<e.first<<" "<<e.second<<endl;
				counting_errors++;
			}
		}
		cout << "The results were checked" << endl;
		cout<< counting_errors<<"	errors"<<endl;
	}
}



void dump_stats() {
	uint64_t total_super_kmers(0);
	uint64_t total_kmers(0);
	uint64_t non_null_buckets(0);
	uint64_t null_buckets(0);
	uint64_t largest_bucket(0);

	//FOREACH BUCKETS
	for (uint64_t i(0); i < bucket_number.value(); ++i) {
		if (bucket_menus[i].size() != 0) {
			for (auto e:bucket_menus[i]) {
				largest_bucket = max(largest_bucket, (uint64_t)bucket_menus[i].size());
				non_null_buckets++;
				total_super_kmers +=bucket_menus[i].size();
				total_kmers += bucket_menus[i].number_kmer();
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




bool count_abundance(const string& line) {
	if (line.size() < k) {
		return false;
	}
	// Init Sequences
	uint64_t min_seq  = (str2num(line.substr(0, minimizer_size))), min_rcseq(rcbc(min_seq, minimizer_size)), min_canon(min(min_seq, min_rcseq));
	// #pragma omp critical(abundance_mini)
	// {
		abundance_mini[min_canon]++;
		if(abundance_mini[min_canon]>65000){
			return true;
		}
		// cout<<min_canon<<" "<<abundance_mini[min_canon]<<endl;
	// }

	uint64_t line_size = line.size();
	for (uint64_t i = 0; i + minimizer_size < line_size; ++i) {
		// UpdateMINIMIZER candidate with the new letter
		updateM(min_seq, line[i + minimizer_size]);
		updateRCM(min_rcseq, line[i + minimizer_size]);
		min_canon = (min(min_seq, min_rcseq));
		// #pragma omp critical(abundance_mini)
		// {
			abundance_mini[min_canon]++;
			if(abundance_mini[min_canon]>65000){
				return true;
			}
			// cout<<min_canon<<" "<<abundance_mini[min_canon]<<endl;
		// }
	}
	return false;
}




void read_fasta_file_ab(const string& filename) {
	zstr::ifstream in(filename);
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
			if(count_abundance(line)){
				break;
			}
			line_count++;
		}
	}
}




void count_line(const string& line) {
	if (line.size() < k) {
		return;
	}
	vector<kmer_full> kmers;
	// Init Sequences
	uint64_t kmer_seq = (str2num(line.substr(0, k))), kmer_rc_seq(rcb(kmer_seq));
	uint64_t min_seq  = (str2num(line.substr(k - minimizer_size, minimizer_size))), min_rcseq(rcbc(min_seq, minimizer_size)), min_canon(min(min_seq, min_rcseq));
	// Init MINIMIZER
	int8_t relative_min_position;
	int64_t minimizer = get_minimizer(kmer_seq, relative_min_position);
	bool multiple_min  = relative_min_position < 0;

	uint8_t position_minimizer_in_kmer;
	if (multiple_min){
		position_minimizer_in_kmer = (uint8_t)(-relative_min_position - 1);
	}else{
		position_minimizer_in_kmer = relative_min_position;
	}

	uint64_t hash_mini = hash64shift(abs(minimizer));
	if(minimizer<0){
		kmers.push_back({31-relative_min_position-minimizer_size, kmer_rc_seq});
	}else{
		kmers.push_back({relative_min_position, kmer_seq});
	}
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
			if(minimizer<0){
				reverse(kmers.begin(),kmers.end());
				minimizer*=-1;
			}
			uint64_t bucketindice = minimizer % bucket_number;
			uint64_t fingerprint=minimizer/ bucket_number;
			omp_set_lock(&MutexWall[bucketindice % 4096]);
			if(buckets_index[bucketindice].count(fingerprint)==0){
				bucket_menus[bucketindice].push_back(Bucket());
				buckets_index[bucketindice][fingerprint]=bucket_menus[bucketindice].size()-1;
			}
			bucket_menus[bucketindice][buckets_index[bucketindice][fingerprint]].add_kmers(kmers);
			omp_unset_lock(&MutexWall[bucketindice % 4096]);
			minimizer                  = (min_canon);
			if(min_canon!=min_seq){minimizer*=-1;}
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
			if(minimizer<0){
				reverse(kmers.begin(),kmers.end());
				minimizer*=-1;
			}
			uint64_t bucketindice = minimizer % bucket_number;
			uint64_t fingerprint=minimizer/ bucket_number;
			omp_set_lock(&MutexWall[bucketindice % 4096]);
			if(buckets_index[bucketindice].count(fingerprint)==0){
				bucket_menus[bucketindice].push_back(Bucket());
				buckets_index[bucketindice][fingerprint]=bucket_menus[bucketindice].size()-1;
			}
			bucket_menus[bucketindice][buckets_index[bucketindice][fingerprint]].add_kmers(kmers);
			omp_unset_lock(&MutexWall[bucketindice % 4096]);
			// Search for the new MINIMIZER in the whole kmer
			minimizer    = get_minimizer(kmer_seq, relative_min_position);
			multiple_min = relative_min_position < 0;
			if (multiple_min){
				position_minimizer_in_kmer = (uint8_t)(-relative_min_position - 1);
			}else{
				position_minimizer_in_kmer = relative_min_position;
			}
			hash_mini = hash64shift(abs(minimizer));
		} else {

			position_minimizer_in_kmer++;
			if (multiple_min){
				relative_min_position--;
			}else{
				relative_min_position++;
			}
		}
		// Normal add of the kmer into kmer list
		if(minimizer<0){
			kmers.push_back({31-relative_min_position-minimizer_size, kmer_rc_seq});
		}else{
			kmers.push_back({relative_min_position, kmer_seq});
		}
	}
	if(minimizer<0){
		reverse(kmers.begin(),kmers.end());
		minimizer*=-1;
	}
	uint64_t bucketindice = minimizer % bucket_number;
	uint64_t fingerprint=minimizer/ bucket_number;
	omp_set_lock(&MutexWall[bucketindice % 4096]);
	if(buckets_index[bucketindice].count(fingerprint)==0){
		bucket_menus[bucketindice].push_back(Bucket());
		buckets_index[bucketindice][fingerprint]=bucket_menus[bucketindice].size()-1;
	}
	bucket_menus[bucketindice][buckets_index[bucketindice][fingerprint]].add_kmers(kmers);
	omp_unset_lock(&MutexWall[bucketindice % 4096]);
}



void read_fasta_file(const string& filename) {
	zstr::ifstream in(filename);
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
			count_line(line);
			line_count++;
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
	int mode(1);
	if (argc > 2) {
		mode = stoi(argv[2]);
	}

	if (mode > 1) {
		check = true;
		cout << "LETS CHECK THE RESULTS" << endl;
	}


	cout << "\n\n\nI count " << argv[1] << endl;
	cout << "Minimizer size:	" << minimizer_size << endl;
	cout << "Number of bucket:	" << bucket_number.value() << endl;
	for (uint64_t i(0); i < 4096; ++i) {
		omp_init_lock(&MutexWall[i]);
	}
	auto start = std::chrono::system_clock::now();
	// read_fasta_file_ab(argv[1]);
	auto end = std::chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;

    cout << "Minimize weight computed elapsed time: " << elapsed_seconds.count() << "s\n";
	start = std::chrono::system_clock::now();
	read_fasta_file(argv[1]);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	cout << "Kmer counted elapsed time: " << elapsed_seconds.count() << "s\n";
	cout << endl;
	if (mode % 2 == 0) {
		dump_counting();
	}
	dump_stats();
	exit(0);
}
