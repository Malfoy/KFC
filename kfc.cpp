#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <unordered_map>
#include "index_min.h"
#include "SolidSampler.hpp"

using namespace std;
using namespace chrono;

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

void clean(string& str) {
	for (unsigned i(0); i < str.size(); ++i) {
		switch (str[i]) {
			case 'a': break;
			case 'A': break;
			case 'c': break;
			case 'C': break;
			case 'g': break;
			case 'G': break;
			case 't': break;
			case 'T': break;
			default: str = ""; return;
		}
	}
	transform(str.begin(), str.end(), str.begin(), ::toupper);
}

#define hash_letter(letter) ((letter >> 1) & 3)

void insert_sequence(SolidSampler& sampler, const string& seq) {
	uint64_t hash = 0;
	uint64_t canon_hash = 0;
	for (unsigned i = 0; i < 31; i++)
		hash = hash << 2 | hash_letter(seq[i]);

	for (unsigned idx = 31; idx < seq.size() /**/; idx++) {
		hash = hash << 2 | hash_letter(seq[idx]);
		canon_hash = min(hash, rcb(hash, 31));
		sampler.insert(reinterpret_cast<uint8_t*>(&canon_hash), sizeof(hash));
	}
}

int main(int argc, char** argv) {
	uint64_t size = (uint64_t(1) << 23);
	if (argc < 2) {
		cerr << "[Fasta file]" << endl;
		exit(0);
	}

	SolidSampler sampler(size);
	string input(argv[1]);

	srand(static_cast<unsigned>(time(NULL)));
	string header, sequence, line;
	ifstream in(input);

	// WE TRY TO FIND THE ABUNDANT KMERS

	high_resolution_clock::time_point t1 = high_resolution_clock::now();

	unsigned nb_sequence = 0;
	while (not in.eof()) {
		if (nb_sequence++ >= 100000) break;
		getline(in, header);
		if (header[0] != '>') {
			continue;
		}
		char c = static_cast<char>(in.peek());
		while (c != '>' and c != EOF) {
			getline(in, line);
			sequence += line;
			c = static_cast<char>(in.peek());
		}
		insert_sequence(sampler, sequence);
		sequence = "";
	}

	vector<uint64_t> abundant_kmer = std::move(sampler.get_kmers());
	sampler.clean();

	// SAMPLING DONE NOW WE DO THE ***EASY*** JOB

	in.clear();
	in.seekg(0, ios::beg);
	cerr << "SAMPLING DONE" << endl;
	cerr << abundant_kmer.size() << " abundant kmer" << endl;
	high_resolution_clock::time_point t2 = high_resolution_clock::now();

	duration<double> time_span1 = duration_cast<duration<double>>(t2 - t1);
	cerr << "It took me " << time_span1.count() << " seconds.\n" << endl;

	cerr << "I BUILD THE INDEX" << endl;
	index_full index(abundant_kmer);
	cerr << "DONE	" << endl;
	high_resolution_clock::time_point t3 = high_resolution_clock::now();

	duration<double> time_span2 = duration_cast<duration<double>>(t3 - t2);
	cerr << "It took me " << time_span2.count() << " seconds.\n" << endl;

	while (not in.eof()) {
		getline(in, header);
		if (header[0] != '>') {
			continue;
		}
		char c = static_cast<char>(in.peek());
		while (c != '>' and c != EOF) {
			getline(in, line);
			sequence += line;
			c = static_cast<char>(in.peek());
		}
		index.insert_seq(sequence);
		//~ cin.get();
		index.clear();
		sequence = "";
	}

	high_resolution_clock::time_point t4 = high_resolution_clock::now();

	duration<double> time_span3 = duration_cast<duration<double>>(t4 - t3);

	cerr << "I FINISHED COUNTING !" << endl;
	cerr << "It took me " << time_span3.count() << " seconds.\n" << endl;

	duration<double> time_span4 = duration_cast<duration<double>>(t4 - t1);
	cerr << "Whole process took " << time_span4.count() << " seconds.\n" << endl;

	cin.get();
	// COUNTING WAS DONE IN RAM I OUTPUT THE RESULT BECAUSE OF THE AMAZING AND POWERFULL SAMPLER
	index.dump_counting(std::cout);
	index.clear(true);
	// MY JOB HERE IS DONE *fly away*
}
