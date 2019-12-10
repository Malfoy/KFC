#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <unordered_map>
#include "SolidSampler.hpp"
#include <gatbl/fastx2kmer.hpp>
#include <gatbl/kmer.hpp>
#include <robin_hood.h>

template<typename Window = gatbl::kmer_window<kmer_t>>
class sampling_parser : public gatbl::details::sequence2kmers_base<sampling_parser<Window>, Window> {
  public:
	using base = gatbl::details::sequence2kmers_base<sampling_parser<Window>, Window>;
	template<typename W>
	sampling_parser(W&& window, size_t memory_size)
	  : base(std::forward<W>(window))
	  , _sampler(memory_size) {}

	std::vector<uint64_t>&& get_aboundant_kmers() && { return std::move(_sampler.get_kmers()); }

  protected:
	friend base;

	void on_kmer() {
		assume(not this->done, "asking for more ?");
		this->done = _sampler.insert(this->get_window().canon());
		if (unlikely(this->done)) {
			std::cerr << "Sampling stopped after " << nb_reads << " reads and " << _sampler.get_nb_kmer_seen() << " kmers seen" << std::endl;
			std::cerr << _sampler;
		}
	}
	void on_run_end() {}
	bool on_chrom() {
		nb_reads++;
		return true;
	}

	SolidSampler _sampler;
	size_t nb_reads = 0;
};

template<typename Window = gatbl::kmer_window<kmer_t>>
class counting_parser : public gatbl::details::sequence2kmers_base<counting_parser<Window>, Window> {
  public:
	using counter_t = uint8_t;

	using base = gatbl::details::sequence2kmers_base<counting_parser<Window>, Window>;
	template<typename W>
	counting_parser(W&& window, std::vector<kmer_t>&& kmers)
	  : base(std::forward<W>(window)) {
		for (kmer_t kmer : kmers) {
			_map.insert({kmer, 0});
		}
		kmers.clear();
	}

	friend std::ostream& operator<<(std::ostream& out, counting_parser& that) {
		using sized_kmer_t = typename Window::sized_kmer_t;
		const gatbl::ksize_t k = that.get_window().size();
		for (auto pair : that._map) {
			out << sized_kmer_t{pair.first, k} << ' ' << unsigned(pair.second) << std::endl;
		}
		return out;
	}

  protected:
	friend base;

	void on_kmer() {
		kmer_t kmer = this->get_window().canon();
		auto it = _map.find(kmer);
		if (it != _map.end()) {
			if (it->second != std::numeric_limits<counter_t>::max()) it->second++;
		}
	}
	void on_run_end() {}
	bool on_chrom() {
		nb_reads++;
		return true;
	}

	robin_hood::unordered_flat_map<kmer_t, uint8_t, gatbl::ReversibleHash> _map;
	size_t nb_reads = 0;
};

int main(int argc, char** argv) {
	if (argc < 3) {
		std::cerr << argv[0] << " input.fa/fq memory_in_MB" << std::endl;
		std::exit(1);
	}

	const std::string input_path = argv[1];
	uint64_t size = size_t(std::stod(argv[2]) * double(size_t(1) << 20));
	const gatbl::ksize_t k = 31;

	// WE TRY TO FIND THE ABUNDANT KMERS

	using clock = std::chrono::high_resolution_clock;
	using std::chrono::duration;
	using std::chrono::duration_cast;

	auto t1 = clock::now();
	std::vector<uint64_t> abundant_kmer;
	{
		sampling_parser<> sampler(k, size);
		sampler.read_fastx(input_path);
		abundant_kmer = std::move(sampler).get_aboundant_kmers();
	}
	auto t2 = clock::now();
	auto time_span1 = duration_cast<duration<double>>(t2 - t1);

	std::cerr << "SAMPLING DONE" << std::endl;
	std::cerr << "It took me " << time_span1.count() << " seconds.\n" << std::endl;

	counting_parser<> counter(k, std::move(abundant_kmer));
	counter.read_fastx(input_path);

	auto t3 = clock::now();
	auto time_span2 = duration_cast<duration<double>>(t3 - t2);
	std::cerr << "It took me " << time_span2.count() << " seconds.\n" << std::endl;

	std::cout << counter;

	// MY JOB HERE IS DONE *fly away*
}
