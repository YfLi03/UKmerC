#include "kmerops.hpp"
#include "logger.hpp"
#include "dnaseq.hpp"
#include "timer.hpp"
#include "paradissort.hpp"
#include <cstring>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <omp.h>

void print_kmer_histogram(const KmerList& kmerlist) {
    #if LOG_LEVEL >= 2

    Logger logger;
    int maxcount = std::accumulate(kmerlist.cbegin(), kmerlist.cend(), 0, [](int cur, const auto& entry) { return std::max(cur, std::get<1>(entry)); });

    std::vector<int> histo(maxcount+1, 0);

    for(size_t i = 0; i < kmerlist.size(); ++i)
    {
        int cnt = std::get<1>(kmerlist[i]);
        assert(cnt >= 1);
        histo[cnt]++;
    }

    std::cout << "#count\tnumkmers" << std::endl;

    for (int i = 1; i < histo.size(); ++i)
    {
        if (histo[i] > 0)
        {
            std::cout << i << "\t" << histo[i] << std::endl;
        }
    }
    std::cout << std::endl;

    #endif
}

int GetKmerOwner(const TKmer& kmer, int ntasks) {
    uint64_t myhash = kmer.GetHash();
    double range = static_cast<double>(myhash) * static_cast<double>(ntasks);
    size_t owner = range / std::numeric_limits<uint64_t>::max();
    assert(owner >= 0 && owner < static_cast<int>(ntasks));
    return static_cast<int>(owner);
}

std::unique_ptr<KmerSeedBuckets>
exchange_kmer(const DnaBuffer& myreads,
     int thr_per_task,
     int max_thr_membounded)
{

    #if LOG_LEVEL >= 3
    Timer timer;
    #endif

    Logger logger;
    std::ostringstream rootlog;

    /* parallel settings */
    omp_set_nested(1);
    int ntasks = omp_get_max_threads() / thr_per_task;
    if (ntasks < 1) {
        ntasks = 1;
        thr_per_task = omp_get_max_threads();
    }
    /* for memory bounded operations in this stage, we use another number of threads */
    int nthr_membounded = std::min(omp_get_max_threads() , max_thr_membounded);

    size_t numreads = myreads.size();     /* Number of locally stored reads */

    #if LOG_LEVEL >= 2
    logger() << ntasks << " \t (thread per task: " << thr_per_task << ")";
    logger.flush("Task num:");
    logger() << nthr_membounded ;
    logger.flush("Thread count used for memory bounded operations:");
    #endif

    #if LOG_LEVEL >= 3
    timer.start();
    #endif

    /* prepare the vars for each task */
    KmerSeedVecs kmerseeds_vecs;    
    for (int i = 0; i < nthr_membounded; i++ ) {
        kmerseeds_vecs.push_back(std::vector<std::vector<KmerSeedStruct>>(ntasks));
    }

    std::vector<KmerParserHandler> parser_vecs;
    for (int i = 0; i < nthr_membounded; i++ ) {
        /* This is a little bit tricky, as we need to avoid the kmerseeds_vecs from reallocating */
        parser_vecs.push_back(KmerParserHandler(kmerseeds_vecs[i]));
    }

    #pragma omp parallel for num_threads(nthr_membounded)
    for (int i = 0; i < nthr_membounded; i++) {
        for(int j = 0; j < ntasks; j++) {
            kmerseeds_vecs[i][j].reserve(size_t(1.1 * numreads / ntasks / nthr_membounded));
        }
    }

    ForeachKmerParallel(myreads, parser_vecs, nthr_membounded);

    #if LOG_LEVEL >= 3
    timer.stop_and_log("K-mer partitioning");
    timer.start();
    #endif

    KmerSeedBuckets* recv_kmerseeds = new KmerSeedBuckets(ntasks);

    #if LOG_LEVEL >= 3
    timer.start();
    #endif

    #pragma omp parallel for num_threads(nthr_membounded)
    for (int j = 0; j < ntasks; j++) {
        for (int i = 0; i < nthr_membounded; i++) {
            (*recv_kmerseeds)[j].insert((*recv_kmerseeds)[j].end(), kmerseeds_vecs[i][j].begin(), kmerseeds_vecs[i][j].end());
        }
    }
    

    #if LOG_LEVEL >= 3
    timer.stop_and_log("Local K-mer format conversion");
    #endif

    /* for 1 process, no MPI communication is necessary*/
    return std::unique_ptr<KmerSeedBuckets>(recv_kmerseeds);
}


std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBuckets>& recv_kmerseeds, int thr_per_task )
{
    Timer timer;

    Logger logger;
    std::ostringstream rootlog;

    int ntasks = omp_get_max_threads() / thr_per_task;
    if (ntasks < 1){
        ntasks = 1;
        thr_per_task = omp_get_max_threads();
    }

    assert(ntasks == recv_kmerseeds->size());

    omp_set_nested(1);
    omp_set_num_threads(ntasks);

    uint64_t task_seedcnt[ntasks];
    uint64_t valid_kmer[ntasks];
    KmerList kmerlists[ntasks];

    #if LOG_LEVEL >= 2
    logger() <<"task num: "<< ntasks << " \tthread per task: " << thr_per_task <<" \ttask size: ";
    for (int i = 0; i < ntasks; i++) {
        task_seedcnt[i] = (*recv_kmerseeds)[i].size();
        logger() << task_seedcnt[i] << " ";
    }
    logger.flush("Parallel tasks for kmer sorting and counting:");
    #endif

    #if LOG_LEVEL >= 3
    timer.start();
    #endif

    // yfli: maybe ask omp to scatter the threads on different places
    #pragma omp parallel 
    {
        int tid = omp_get_thread_num();
        paradis::sort<KmerSeedStruct, TKmer::NBYTES>((*recv_kmerseeds)[tid].data(), (*recv_kmerseeds)[tid].data() + task_seedcnt[tid], thr_per_task);
    

    #if LOG_LEVEL >= 3
    }
    /* this implicit barrier can be eliminated when not debugging */
    timer.stop_and_log("Shared memory parallel K-mer sorting");
    timer.start();
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
    #endif

        kmerlists[tid].reserve(uint64_t(task_seedcnt[tid] / LOWER_KMER_FREQ));    // This should be enough
        TKmer last_mer = (*recv_kmerseeds)[tid][0].kmer;
        uint64_t cur_kmer_cnt = 1;
        valid_kmer[tid] = 0;

        for(size_t idx = 1; idx < task_seedcnt[tid]; idx++) {
            TKmer cur_mer = (*recv_kmerseeds)[tid][idx].kmer;
            if (cur_mer == last_mer) {
                cur_kmer_cnt++;
            } else {
                if (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ) {

                    kmerlists[tid].push_back(KmerListEntry());
                    KmerListEntry& entry    = kmerlists[tid].back();

                    TKmer& kmer             = std::get<0>(entry);
                    int& count              = std::get<1>(entry);

                    count = cur_kmer_cnt;
                    kmer = last_mer;

                    valid_kmer[tid]++;
                }

                cur_kmer_cnt = 1;
                last_mer = cur_mer;
            }
        }

        /* deal with the last kmer */
        if (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ) {
            kmerlists[tid].push_back(KmerListEntry());
            KmerListEntry& entry         = kmerlists[tid].back();

            TKmer& kmer             = std::get<0>(entry);
            int& count              = std::get<1>(entry);

            count = cur_kmer_cnt;
            kmer = last_mer;
            valid_kmer[tid]++;
        }
    }

    #if LOG_LEVEL >= 3
    timer.stop_and_log("K-mer counting");
    timer.start();
    #endif

    #if LOG_LEVEL >= 2
    for (int i = 0; i < ntasks; i++) {
        logger() << valid_kmer[i] << " ";
    }
    logger.flush("Valid kmer from tasks:");
    #endif


    uint64_t valid_kmer_total = std::accumulate(valid_kmer, valid_kmer + ntasks, (uint64_t)0);
    KmerList* kmerlist = new KmerList();
    kmerlist->reserve(valid_kmer_total);

    // yfli: actually we don't need to copy it into one single kmerlist
    for (int i = 0; i < ntasks; i++) {
        kmerlist->insert(kmerlist->end(), kmerlists[i].begin(), kmerlists[i].end());
    }

    logger() << valid_kmer_total;
    logger.flush("Valid kmer for process:");

    #if LOG_LEVEL >= 3
    timer.stop_and_log("K-mer copying");
    #endif

    return std::unique_ptr<KmerList>(kmerlist);
}