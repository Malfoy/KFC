#include <gatb/gatb_core.hpp>
#include <gatb/kmer/impl/PartitionsCommand.hpp>

using namespace std;

class DSK : public gatb::core::tools::misc::impl::Tool
{
public:

    /** Constructor. */
    DSK ();

private:

    /** \copydoc Tool::execute. */
    void  execute ();
};

struct Parameter
{
    Parameter (DSK& dsk, IProperties* props) : dsk(dsk), props(props) {}
    DSK&         dsk;
    IProperties* props;
};

template<size_t span=KMER_DEFAULT_SPAN>
class CountProcessorDumpTxt : public CountProcessorAbstract<span>
{
public:

    /** Shortcuts. */
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::Type Type;

    /** Constructor */
    CountProcessorDumpTxt (
        size_t                                  kmerSize,
        //system::ISynchronizer*                  synchronizer = 0,
        //tools::storage::impl::Partition<Count>* solidCounts  = 0,
        size_t                                  nbPartsPerPass = 0
    )
        : _kmerSize(kmerSize), _nbPartsPerPass(nbPartsPerPass) 
        //_synchronizer(0), 
         // _solidCounts(0), _solidKmers(0)
    {
    }

    /** Destructor */
    virtual ~CountProcessorDumpTxt ()
    {
    }

    /********************************************************************/
    /*   METHODS CALLED ON THE PROTOTYPE INSTANCE (in the main thread). */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::begin */
    void begin (const Configuration& config)
    {
        /** We remember the number of partitions for one pass. */
        _nbPartsPerPass = config._nb_partitions;

        /** We compute the number of partitions. */
        size_t nbTotalPartitions = config._nb_partitions * config._nb_passes;
    }

    /** \copydoc ICountProcessor<span>::clones */
    CountProcessorAbstract<span>* clone ()
    {
        /** Note : we share the synchronizer for all the clones. */
        return new CountProcessorDumpTxt ( _kmerSize, _nbPartsPerPass);
    }

    /** \copydoc ICountProcessor<span>::finishClones */
    void finishClones (std::vector<ICountProcessor<span>*>& clones)
    {
        for (size_t i=0; i<clones.size(); i++)
        {
            /** We have to recover type information. */
            if (CountProcessorDumpTxt* clone = dynamic_cast<CountProcessorDumpTxt*> (clones[i]))
            {
                for (std::map<std::string,size_t>::iterator it = clone->_namesOccur.begin(); it != clone->_namesOccur.end(); ++it)
                {
                    this->_namesOccur[it->first] += it->second;
                }
            }
        }
    }

    /********************************************************************/
    /*   METHODS CALLED ON ONE CLONED INSTANCE (in a separate thread).  */
    /********************************************************************/

    /** \copydoc ICountProcessor<span>::beginPart */
    void beginPart (size_t passId, size_t partId, size_t cacheSize, const char* name)
    {
        /** We get the actual partition idx in function of the current partition AND pass identifiers. */
        size_t actualPartId = partId + (passId * _nbPartsPerPass);

        /** We get a handle on the current solid bag (with a cache). */
        //setSolidKmers (new tools::collections::impl::BagCache<Count> (& (*_solidCounts)[actualPartId], cacheSize, _synchronizer));

        /** We update some stats (want to know how many "hash" or "vector" partitions we use). */
        _namesOccur[name] ++;
    }

    /** \copydoc ICountProcessor<span>::endPart */
    void endPart (size_t passId, size_t partId)
    {
        /** We flush the current collection for the partition just processed. */
        //_solidKmers->flush();
    }

    /** \copydoc ICountProcessor<span>::process */
    bool process (size_t partId, const Type& kmer, const CountVector& count, CountNumber sum)
    {
        //this->_solidKmers->insert (Count(kmer,sum));
        
        // TODO output the kmer here
        std::cout << kmer << sum << std::endl;
        return true;
    }

private:

    size_t _kmerSize;

    size_t _nbPartsPerPass;

    std::map<std::string,size_t> _namesOccur;
};



template<size_t span> struct Functor  {  void operator ()  (Parameter parameter)
{
    DSK&         dsk   = parameter.dsk;
    IProperties* props = parameter.props;

    /** We get a handle on tha bank. */
    IBank* bank = Bank::open(props->getStr("-file"));
    LOCAL (bank);


    // TODO change all of this 
    // ---------------snip------------------
	uint nbCores = 8;
	uint nbCores_per_partition = 1;
	uint kmerSize = 31;
    uint max_memory = 4000;
    
    uint current_core = 0;

    uint nb_partitions = 100;
    uint minim_size = 10;
    PartiInfo<5> pInfo (nb_partitions, minim_size);// not sure if useful

    // FIXME put the number of items in the bucket here
    vector<size_t> nbItemsPerBankPerPart;
    for (size_t i=0; i<nb_partitions; i++)
    {
        nbItemsPerBankPerPart.push_back (1000000);
    }

    u_int64_t mem = (max_memory*MBYTE)/nbCores;
    typedef typename Kmer<span>::Count  Count;
	size_t cacheSize = std::min ((u_int64_t)(200*1000), mem/(50*sizeof(Count)));

    std::string _tmpStorageName_superK; // TODO maybe needs to be set
    gatb::core::tools::storage::impl::SuperKmerBinFiles* _superKstorage = new SuperKmerBinFiles(_tmpStorageName_superK,"superKparts", nb_partitions) ;

    // ---------------snip------------------
    
	MemAllocator pool (nbCores);
    
	uint pass = 0; // probably shouldn't change
    
    gatb::core::tools::misc::impl::TimeInfo _fillTimeInfo;

    // si on voulait utiliser creer un histo, suffit de decommenter:
    // (et c'est aussi combinable avec le dump des counts)
    /* gatb::core::kmer::impl::CountProcessorHistogram<span>* hist_processor =  new CountProcessorHistogram<span> (
            0,
            1000, //params->getInt(STR_HISTOGRAM_MAX),
            1, // params->getInt(STR_KMER_ABUNDANCE_MIN_THRESHOLD),
			false, // using_histo_2D,
			true, // using_histo_1D,
			"histo2D", //histo2Dstorage_filename,
			"histo1D" //histo1Dstorage_filename
        );
    */
    /*
    typedef ICountProcessor<span> CountProcessor;
    vector<CountProcessor*> clones;
	CountProcessor* processorClone = processor->clone ();
	processorClone->use();
	clones.push_back (processorClone);
    */
    CountProcessorDumpTxt<span>* dump_processor = new CountProcessorDumpTxt<span>(kmerSize, nb_partitions);
    
    gatb::core::tools::dp::IteratorListener* _progress;
    _progress->init ();

    ICommand* cmd = new PartitionsByVectorCommand<span> (dump_processor, cacheSize, _progress, _fillTimeInfo,
											   pInfo, pass, current_core, nbCores_per_partition, kmerSize, pool, nbItemsPerBankPerPart, _superKstorage);

	cmd->execute();

} };

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
DSK::DSK () : Tool ("dsk")
{
    /** We add options specific to DSK (most important at the end). */
    getParser()->push_back (SortingCountAlgorithm<>::getOptionsParser(), 1);

    /** We rename the input option. */
    if (IOptionsParser* input = getParser()->getParser (STR_URI_INPUT))  {  input->setName (STR_URI_FILE);  }
}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
void DSK::execute ()
{
    /** we get the kmer size chosen by the end user. */
    size_t kmerSize = getInput()->getInt (STR_KMER_SIZE);

    /** We launch dsk with the correct Integer implementation according to the choosen kmer size. */
    Integer::apply<Functor,Parameter> (kmerSize, Parameter (*this, getInput()));
}

/********************************************************************************/

int main (int argc, char* argv[])
{
    // We define a try/catch block in case some method fails (bad filename for instance)
    try
    {
        /** We execute dsk. */
        DSK().run (argc, argv);
    }

    catch (OptionFailure& e)
    {
        return e.displayErrors (std::cout);
    }

    catch (Exception& e)
    {
        cerr << "EXCEPTION: " << e.getMessage() << endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
