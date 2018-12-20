#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include "index_min.h"

#include "BloomFilter.hpp"




using namespace std;


string getLineFasta(ifstream* in){
	string line,result;
	getline(*in,line);
	char c=in->peek();
	while(c!='>' and c!=EOF){
		getline(*in,line);
		result+=line;
		c=in->peek();
	}
	return result;
}


void clean(string& str){
	for(uint i(0); i< str.size(); ++i){
		switch(str[i]){
			case 'a':break;
			case 'A':break;
			case 'c':break;
			case 'C':break;
			case 'g':break;
			case 'G':break;
			case 't':break;
			case 'T':break;
			default: str="";return;
		}
	}
	transform(str.begin(), str.end(), str.begin(), ::toupper);
}

#define hash_letter(letter) ((letter >> 1) & 0b11)

void insert_sequence(BloomFilter& bf, const string& seq){
	uint64_t hash = 0;
	for (uint i=0 ; i<31 ; i++)
		hash = hash << 2 | hash_letter(seq[i]);

	for (uint idx=31 ; idx<31+20/*seq.size()/**/ ; idx++) {
		hash = hash << 2 | hash_letter(seq[idx]);
		bf.add((uint8_t *)&hash, sizeof(hash));
	}
}



int main(int argc, char ** argv){
	uint64_t size = 50;
	// size <<= 30;
	BloomFilter bf(size, 3);

	if(argc<2){
		cout<<"[Fasta file]"<<endl;
		exit(0);
	}
	string input(argv[1]);

	srand (time(NULL));
	string header, sequence,line;
	ifstream in(input);


	while(not in.eof()){
		getline(in,header);
		if(header[0]!='>'){continue;}
		char c=in.peek();
		while(c!='>' and c!=EOF){
			getline(in,line);
			sequence+=line;
			c=in.peek();
		}
		insert_sequence(bf, sequence);
		cout << bf;
		// TODO: remove this
		return 0;
		sequence="";
	}
	vector<uint64_t> abundant_kmer;
	index_full index(abundant_kmer);
}
