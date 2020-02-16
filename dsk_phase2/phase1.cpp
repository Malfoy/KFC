#include <gatb/gatb_core.hpp>
#include <gatb/kmer/impl/RepartitionAlgorithm.hpp>
#include <gatb/kmer/impl/Configuration.hpp>
#include <gatb/kmer/impl/ConfigurationAlgorithm.hpp>
#include <gatb/kmer/impl/SortingCountAlgorithm.cpp> // for FillPartitions
#include <gatb/kmer/impl/Model.hpp>

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



template<size_t span> struct Functor  {  void operator ()  (Parameter parameter)
{
    DSK&         dsk   = parameter.dsk;
    IProperties* props = parameter.props;

    /** We get a handle on tha bank. */
    IBank* bank = Bank::open(props->getStr("-file"));
    LOCAL (bank);


    // TODO change all of this 
    // ---------------snip------------------
	uint nbCores = 1;//8;
	uint nbCores_per_partition = 1;
	uint kmerSize = 31;
    uint max_memory = 4000;
    uint minim_size = 10;

    /** We update the message of the progress bar. */
    // if we have an estimation of file volume, put it here for the progress bar
    //
    /*int nb_iterations= (1 + nbCores) * _config._volume * MBYTE / sizeof(Type);
    gatb::core::tools::dp::IteratorListener* _progress(new ProgressSynchro (
                                            new Progress (nb_iterations, "format0"),
                                                                            System::thread().newSynchronizer())
                                );
    */
    // quiet progress
    gatb::core::tools::dp::IteratorListener* _progress(new ProgressSynchro (
                                            new IteratorListener (),
                                           System::thread().newSynchronizer()));
 
    _progress->init();
    _progress->setMessage ("DSK phase 1");

    /** We create a kmer model; using the frequency order if we're in that mode */
    uint32_t* freq_order = NULL;
    
    // many shenanigans to obtain a repartition of partitions
    // could perhaps use antoine's minimizer repartition trick instead, but we didn't know it when making gatb-core
    gatb::core::kmer::impl::ConfigurationAlgorithm<span> configAlgo (bank, props);
    configAlgo.execute();
    gatb::core::kmer::impl::Configuration _config = configAlgo.getConfiguration();

    string prefix = basename(props->getStr("-file").c_str());
    gatb::core::tools::storage::impl::Storage* storage = StorageFactory(gatb::core::tools::storage::impl::STORAGE_FILE).create (prefix+"_storage", true, false);
    RepartitorAlgorithm<span> repart (
                bank, 
                storage->getGroup("minimizers"), 
                _config,
                nbCores
                );
        repart.execute ();
    Repartitor*             repartitor = new Repartitor(storage->getGroup("minimizers"));

    string _tmpStorageName_superK = prefix+"_superK_partitions";
    std::cout << "nb_partitions determined by ConfigurationAlgorithm: " << _config._nb_partitions << std::endl;;
    std::cout << "nb cores from ConfigurationAlgorithm: " <<  _config._nbCores << " " << _config._nbCores_per_partition << " " << _config._nb_partitions_in_parallel <<std::endl;
    uint nb_partitions = _config._nb_partitions;

    gatb::core::tools::storage::impl::SuperKmerBinFiles*  _superKstorage = new SuperKmerBinFiles(_tmpStorageName_superK,"superKparts", nb_partitions) ;

    /** We may have to retrieve the minimizers frequencies computed in the RepartitorAlgorithm. */
    int minimizerType = 1;
    if (minimizerType == 1)  {  freq_order = repartitor->getMinimizerFrequencies ();  }

    typedef typename Kmer<span>::Type                                       Type;
    typedef typename Kmer<span>::ModelCanonical                         ModelCanonical;
    typedef typename Kmer<span>::template ModelMinimizer <ModelCanonical>   Model;
    Model model( kmerSize, minim_size, typename gatb::core::kmer::impl::Kmer<span>::ComparatorMinimizerFrequencyOrLex(), freq_order);

    /** We have to reinit the progress instance since it may have been used by SampleRepart before. */
    _progress->init();


    /** We create the sequences iterator. */
    Iterator<Sequence>* itSeq = bank->iterator();
    LOCAL (itSeq);


    BankStats _bankStats;
    PartiInfo<5> pInfo (nb_partitions, minim_size);
    
    { // these brackets are needed so that FillPartitions gets deleted
     // otherwise superKstorage->closefiles() gets called before writing gets to be done
        auto fillpartitions = FillPartitions<span,true>(
                model, 1 /* nb_passes */, 0 /* pass */, nb_partitions,
                _config._nb_cached_items_per_core_per_part, _progress, _bankStats,
                nullptr /*_tmpPartitions*/ /* unneeded with new FillPartition with superkmers*/
                , *repartitor, pInfo, _superKstorage);

        /* We fill the partitions. Each thread will read synchronously and will
         * call FillPartitions in a synchronous way (in order to have global
         * BanksStats correctly computed). *//*
                                                size_t groupSize = 1000;
                                                bool deleteSynchro = true;
                                                getDispatcher()->iterate(
                                                itSeq,
                                                fillpartitions,
                                                groupSize, deleteSynchro);
                                                */

        for (itSeq->first(); !itSeq->isDone(); itSeq->next())
        {
            fillpartitions(itSeq->item());
        }

        // GR: close the input bank here with call to finalize
        itSeq->finalize();
    }

    _superKstorage->flushFiles();
    _superKstorage->closeFiles();

    _superKstorage->saveInfoFile(_tmpStorageName_superK);
    pInfo.saveInfoFile(_tmpStorageName_superK);
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
