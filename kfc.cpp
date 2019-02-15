#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <unordered_map>
#include "index_min.h"
#include "SolidSampler.hpp"



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






void insert_sequence(SolidSampler& sampler, const string& seq){
	uint64_t hash = 0;
	uint64_t canon_hash = 0;
	for (uint i=0 ; i<31 ; i++)
		hash = hash << 2 | hash_letter(seq[i]);

	for (uint idx=31 ; idx<seq.size()/**/ ; idx++) {
		hash = hash << 2 | hash_letter(seq[idx]);
		canon_hash=min(hash,rcb(hash,31));
		sampler.insert((uint8_t *)&canon_hash, sizeof(hash));
	}
}



int main(int argc, char ** argv){
  uint64_t size = ((uint64_t)1 << 33);
	if(argc<2){
		cout<<"[Fasta file]"<<endl;
		exit(0);
	}

	SolidSampler sampler(size);
	string input(argv[1]);

	srand (time(NULL));
	string header, sequence,line;
	ifstream in(input);

	//WE TRY TO FIND THE ABUNDANT KMERS

	uint nb_sequence=0;
	while(not in.eof()){
		if (nb_sequence++ >= 10000)
			break;
		getline(in,header);
		if(header[0]!='>'){continue;}
		char c=in.peek();
		while(c!='>' and c!=EOF){
			getline(in,line);
			sequence+=line;
			c=in.peek();
		}
		insert_sequence(sampler, sequence);
		sequence="";
	}

	vector<uint64_t> abundant_kmer = *(sampler.get_kmers());
	sampler.clean();

	//SAMPLING DONE NOW WE DO THE ***EASY*** JOB

    in.clear();
	in.seekg(0, ios::beg);

	cout<<abundant_kmer.size()<<" abundant kmer"<<endl;
	cout<<"I BUILD THE INDEX"<<endl;
	index_full index(abundant_kmer);
	cout<<"DONE	"<<endl;
	while(not in.eof()){
		getline(in,header);
		if(header[0]!='>'){continue;}
		char c=in.peek();
		while(c!='>' and c!=EOF){
			getline(in,line);
			sequence+=line;
			c=in.peek();
		}
		index.insert_seq(line);
		//~ cin.get();
		index.clear();
		sequence="";
	}
	cout<<"I FINISHED COUNTING !"<<endl;
	cin.get();
	//COUNTING WAS DONE IN RAM I OUTPUT THE RESULT BECAUSE OF THE AMAZING AND POWERFULL SAMPLER
	index.dump_counting();
	index.clear(true);
	//MY JOB HERE IS DONE *fly away*
}









