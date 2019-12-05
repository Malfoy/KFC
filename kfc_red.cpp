#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <unordered_map>

using namespace std;

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


struct SKC{
  string sk;
  vector<uint8_t> counts;
};


uint64_t k=(31);

uint64_t minimizer_size(9);

uint64_t minimizer_number((uint64_t)1<<(2*minimizer_size));



  uint32_t revhash ( uint32_t x ) {
	x = ( ( x >> 16 ) ^ x ) * 0x2c1b3c6d;
	x = ( ( x >> 16 ) ^ x ) * 0x297a2d39;
	x = ( ( x >> 16 ) ^ x );
	return x;
}



  uint32_t unrevhash ( uint32_t x ) {
	x = ( ( x >> 16 ) ^ x ) * 0x0cf0b109; // PowerMod[0x297a2d39, -1, 2^32]
	x = ( ( x >> 16 ) ^ x ) * 0x64ea2d65;
	x = ( ( x >> 16 ) ^ x );
	return x;
}



  uint64_t revhash ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x ) * 0xD6E8FEB86659FD93;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



uint64_t unrevhash ( uint64_t x ) {
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x ) * 0xCFEE444D8B59A89B;
	x = ( ( x >> 32 ) ^ x );
	return x;
}



uint64_t rcbc(uint64_t in, uint64_t n){
  //assume(n <= 32, "n=%u > 32", n);
  // Complement, swap byte order
  uint64_t res = __builtin_bswap64(in ^ 0xaaaaaaaaaaaaaaaa);
  // Swap nuc order in bytes
  const uint64_t c1 = 0x0f0f0f0f0f0f0f0f;
  const uint64_t c2 = 0x3333333333333333;
  res               = ((res & c1) << 4) | ((res & (c1 << 4)) >> 4); // swap 2-nuc order in bytes
  res               = ((res & c2) << 2) | ((res & (c2 << 2)) >> 2); // swap nuc order in 2-nuc

  // Realign to the right
  res >>= 64 - 2 * n;
  return res;
}



uint64_t canonize(uint64_t x,uint64_t n){
	return min(x,rcbc(x,n));
}


uint64_t get_minimizer(uint64_t seq){
	uint64_t mini,mmer,hash_mini(-1);
	mmer=seq%minimizer_number;
	mini=mmer=canonize(mmer,minimizer_size);
	for(uint64_t i(1);i<=k-minimizer_size;i++){
		seq>>=2;
		mmer=seq%minimizer_number;
		mmer=canonize(mmer,minimizer_size);
		uint64_t hash = (unrevhash(mmer));
		if(hash_mini>hash){
			mini=mmer;
			hash_mini=hash;
		}
	}
	return revhash((uint64_t)mini)%minimizer_number;
}


uint64_t nuc2int(char c){
	return (c/2)%4;
}



uint64_t nuc2intrc(char c){
	return ((c/2)%4)^2;
}



uint64_t offsetUpdateAnchor((uint64_t)1<<62);



inline void updateK(uint64_t& min, char nuc){
	min<<=2;
	min+=nuc2int(nuc);
	min%=offsetUpdateAnchor;
}



inline void updateRCK(uint64_t& min, char nuc){
	min>>=2;
	min+=(nuc2intrc(nuc)<<(2*k-2));
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


int64_t kmer_in_super_kmer(const string& super_kmer,const string& kmer){
	size_t found = super_kmer.find(kmer);
    if (found != string::npos){
			return found;
		}
		found = super_kmer.find(revComp(kmer));
	    if (found != string::npos){
				return found;
			}
		return -1;
}


int compact(string& super_kmer,const string& kmer){
	uint64_t ssk(super_kmer.size()),sk(kmer.size());
	if(sk==0){return -1;}
	if(ssk==0){
		super_kmer=kmer;
		return 0;
	}
	string rc_k(revComp(kmer)),end_sk(super_kmer.substr(ssk-k+1,k-1)), beg_k(kmer.substr(0,k-1));
	if(end_sk==beg_k){
		super_kmer=super_kmer+(kmer.substr(k-1));
		return super_kmer.size()-k+1;
	}
	string beg_rck(rc_k.substr(0,k-1));
	if(end_sk==beg_rck){
		super_kmer=super_kmer+(rc_k.substr(k-1));
		return super_kmer.size()-k+1;
	}
	return -1;
}


uint64_t str2num(const string& str){
uint64_t res(0);
for(uint64_t i(0);i<str.size();i++){
	res<<=2;
	res+=(str[i]/2)%4;
}
return res;
}



void dump_count(const SKC& skc){
	// cout<<"DP"<<endl;
	for(uint64_t i(0);i<skc.counts.size();++i){
		// cout<<i<<endl;
		cout<<skc.sk.substr(i,k)<<" "<<(uint64_t)skc.counts[i]<<endl;
	}
}


void insert_kmer(const string& str_kmer, vector<vector<SKC>>& skc){
	// cout<<"INSERT: "<<str_kmer<<endl;
  uint64_t min(get_minimizer(str2num(str_kmer)));
	// cout<<min<<" "<<skc.size()<<endl;
	for(uint64_t i(0); i < skc[min].size();i++){
	 	int64_t pos(kmer_in_super_kmer(skc[min][i].sk,str_kmer));
    if(pos>=0){
      skc[min][i].counts[pos]++;
			return;
    }
  }

  for(uint64_t i(0); i < skc[min].size();i++){
	 	int64_t pos(compact(skc[min][i].sk,str_kmer));
    if(pos>	0){
      if(skc[min][i].counts.size()<pos){
        skc[min][i].counts.resize(pos,0);
      }
      skc[min][i].counts[pos-1]++;
			return;
    }
  }
	skc[min].push_back({str_kmer,{1}});
}




void dump_counting(vector<vector<SKC>>& buckets){
	for(uint64_t i(0);i< buckets.size();++i) {
		for(uint64_t ii(0);ii<buckets[i].size();++ii){
			// cout<<i<<" "<<ii<<" ";
			dump_count(buckets[i][ii]);
		}
	}
}


string getCanonical(const string& str){
 return (min(str,revComp(str)));
}
	uint64_t nb_kmer_read(0);
void count_line(const string& line, vector<vector<SKC>>& buckets){
  string str_kmer;

  for(uint64_t i=0;i+k<=line.size();++i){
    str_kmer=getCanonical(line.substr(i,k));
    insert_kmer(str_kmer,buckets);
		nb_kmer_read++;
  }
	// nb_kmer_read+= line.size()-k+1;
	// cout<<nb_kmer_read<<endl;
	// cin.get();
}



void read_fasta_file(const string& filename,vector<vector<SKC>>& buckets){
  ifstream in(filename);
  string line;
	uint64_t line_count;
  while(in.good()){
    line=getLineFasta(&in);
    count_line(line,buckets);
		line_count++;
		if(line_count%1000==0){
			cout<<"-"<<flush;
		}
  }
	// cout<<nb_kmer_read<<endl;
	// cin.get();
}




int main(int argc, char ** argv){
	if(argc<2){
		cout<<"[fasta file]"<<endl;
		exit(0);
	}
	vector<vector<SKC>> buckets(minimizer_number);
  read_fasta_file(argv[1],buckets);
	dump_counting(buckets);
  return 0;
}
