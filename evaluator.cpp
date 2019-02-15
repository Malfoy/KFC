#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <unordered_set>
#include "sparsepp/spp.h"
#include <omp.h>
#include <sstream>

using spp::sparse_hash_map;



//~ uint64_t xs(uint64_t y){
	//~ y^=(y<<13); y^=(y>>17);y=(y^=(y<<15)); return y;
//~ }

using namespace std;


string intToString(uint64_t n){
	if(n<1000){
		return to_string(n);
	}
	string end(to_string(n%1000));
	if(end.size()==3){
		return intToString(n/1000)+","+end;
	}
	if(end.size()==2){
		return intToString(n/1000)+",0"+end;
	}
	return intToString(n/1000)+",00"+end;
}



char revCompChar(char c) {
	switch (c) {
		case 'A': return 'T';
		case 'C': return 'G';
		case 'G': return 'C';
	}
	return 'A';
}



string revComp(const string& s){
	string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = revCompChar(s[i]);
	}
	return rc;
}



string getCanonical(const string& str){
	return (min(str,revComp(str)));
}



uint64_t str2num(const string& str){
	uint64_t res(0);
	for(uint64_t i(0);i<str.size();i++){
		res<<=2;
		switch (str[i]){
			case 'A':res+=0;break;
			case 'C':res+=1;break;
			case 'G':res+=2;break;
			default:res+=3;break;
		}
	}
	return res;
}

vector<string> split(const string &s, char delim){
    stringstream ss(s);
    string item;
    vector<string> elems;
    while (getline(ss, item, delim)) {
        elems.push_back(move(item)); 
    }
    return elems;
}


int main(int argc, char ** argv){
	if(argc<4){
		cout<<"[kmer count file to evaluate] [reference kmer count file] [k value] [core number] [n for 2^n pass]"<<endl;
		exit(0);
	}
	auto start = chrono::system_clock::now();
	string inputUnitig(argv[1]);
	string inputRef(argv[2]);
	uint k(stoi(argv[3]));
	uint n(0);
	uint nb_cores(stoi(argv[4]));
	if(argc>5){
		 n=(stoi(argv[5]));
	}
	uint nbHash=1<<n;
	cout<<"I will perform "<<nbHash<<" pass"<<endl;
	srand (time(NULL));

	ifstream inRef(inputRef),inUnitigs(inputUnitig);
	if(not inRef.good() or not inUnitigs.good()){
		cout<<"Problem with files opening"<<endl;
		exit(1);
	}
	uint64_t FP(0),TP(0),FN(0),size(0),number(0),genomicKmersNum(0), TPcount(0), undercount(0), overcount(0);
	omp_lock_t lock[1024];
	for (int i=0; i<1024; i++)
        omp_init_lock(&(lock[i]));

	for(uint HASH(0);HASH<nbHash;++HASH){
		vector<sparse_hash_map<string, uint64_t>> genomicKmers;
		genomicKmers.resize(1024);

		#pragma omp parallel num_threads(nb_cores)
		{
			string ref, header,canon;
			uint64_t headerInt;
			vector<string> headerElements;
			while(not inRef.eof()){
				#pragma omp critical(dataupdate)
				{
					getline(inRef,header);
					getline(inRef,ref);
				}
				if(not ref.empty() and not header.empty()){
					for(uint i(0);i+k<=ref.size();++i){
						canon=(getCanonical(ref.substr(i,k)));
						uint64_t num((str2num(canon)));
						if(num%nbHash==HASH){
							headerElements = split(header,'>');
							headerInt = stoi(headerElements[1]);
							uint64_t num2( (num/nbHash)%1024);
							omp_set_lock(&(lock[num2]));
							genomicKmers[num2][canon]=headerInt;
							omp_unset_lock(&(lock[num2]));
							#pragma omp atomic
							genomicKmersNum++;
						}
					}
				}
			}
		}

		#pragma omp parallel num_threads(nb_cores)
		{
			string ref, header,canon;
			uint64_t headerIntResult;
			vector<string> headerElementsResult;
			while(not inUnitigs.eof()){
				#pragma omp critical(dataupdate)
				{
					getline(inUnitigs,header);
					getline(inUnitigs,ref);
				}
				if(not ref.empty() and not header.empty()){
					headerElementsResult = split(header,'>');
					headerIntResult = stoi(headerElementsResult[1]);
					#pragma omp atomic
					size+=ref.size();
					#pragma omp atomic
					number++;
					for(uint i(0);i+k<=ref.size();++i){
						canon=getCanonical(ref.substr(i,k));
						uint64_t num((str2num(canon)));
						if(num%nbHash==HASH){
							if(genomicKmers[(num/nbHash)%1024].count(canon)==0){
								#pragma omp atomic
								FP++;
							}else{
								#pragma omp atomic
								TP++;
								if(genomicKmers[(num/nbHash)%1024][canon] == headerIntResult){
									#pragma omp atomic
									TPcount++;
								} else if (genomicKmers[(num/nbHash)%1024][canon] > headerIntResult){
									#pragma omp atomic
									overcount++;
								} else {
									#pragma omp atomic
									undercount++;
								}
							}
						}
					}
				}
			}
		}
		if(HASH==0){
			//~ cout<<"k-mer number: "<<intToString(number)<< " Total size: "<<intToString(size)<<" Mean: "<<intToString(size/number)<<endl;
			//~ cout<<"Genomic kmer in the reference: "<<intToString(genomicKmersNum)<<endl;
		}
		cout <<"*******************" << intToString(genomicKmersNum)<< " " << TP << endl;
		FN=genomicKmersNum-TP;
		if(HASH!=nbHash-1){
			cout<<"PARTIAL RESULTS:"<<endl;
			cout<<"-------- Kmers ---------------"<<endl;
			cout<<"True positive (kmers in the results and the reference) 		GOOD kmers:	"<<intToString(TP)<<endl;
			cout<<"False positive (kmers in the results and NOT in the reference)	ERRONEOUS kmers:	"<<intToString(FP)<<endl;
			cout<<"False Negative (kmers NOT in the result but in the reference)	MISSING kmers:	"<<intToString(FN)<<endl;
			cout<<"Erroneous kmer rate (*10,000): "<<(double)10000*FP/(FP+TP)<<endl;
			cout<<"Missing kmer rate (*10,000): "<<(double)10000*FN/genomicKmersNum<<endl;
			cout<<"------- Counts --------------- "<<endl;
			cout<< "True positive counts (same counts in results and reference) 		:	"<<intToString(TPcount)<<endl;
			cout<< "k-mers over-estimated counts (kmers with higher counts in results than in ref) 		:	"<<intToString(overcount)<<endl;
			cout<< "k-mers with under-estimated counts (kmers with lower counts in results than in ref) 		:	"<<intToString(undercount)<<endl;
		}
		inUnitigs.clear();
		inUnitigs.seekg(0, std::ios::beg);
		inRef.clear();
		inRef.seekg(0, std::ios::beg);
	}
	cout<<endl<<"FINAL RESULTS:"<<endl;
	//~ cout<<genomicKmersNum<<" "<<TP<<endl;
	FN=genomicKmersNum-TP;

	cout<<"-------- Kmers ---------------"<<endl;
	cout<<"True positive k-mers(k-mers in the results and the reference) :	"<<intToString(TP)<<endl;
	cout<<"False positive k-mers(k-mers in the results and NOT in the reference):	"<<intToString(FP)<<endl;
	cout<<"False Negative k-mers(k-mers NOT in the result but in the reference):	"<<intToString(FN)<<endl;
	cout<<"Erroneous kmer rate (*10,000): "<<(double)10000*FP/(FP+TP)<<endl;
	cout<<"Missing kmer rate (*10,000): "<<(double)10000*FN/genomicKmersNum<<endl;
	cout<<"------- Counts --------------- "<<endl;
	cout<< "True positive counts (same counts in results and reference):	"<<intToString(TPcount)<<endl;
	cout<< "k-mers over-estimated counts (k-mers with higher counts in results than in ref):	"<<intToString(overcount)<<endl;
	cout<< "k-mers with under-estimated counts (k-mers with lower counts in results than in ref):	"<<intToString(undercount)<<endl;

	auto end = chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
    time_t end_time = chrono::system_clock::to_time_t(end);

    cout << "\nFinished computation at " << ctime(&end_time)<< "Elapsed time: " << elapsed_seconds.count() << "s\n";

}
