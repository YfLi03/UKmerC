#ifndef KMEROPS_HPP
#define KMEROPS_HPP

#include "kmer.hpp"
#include "timer.hpp"
#include "dnaseq.hpp"
#include "dnabuffer.hpp"
#include "logger.hpp"
#include <omp.h>


typedef uint32_t PosInRead;
typedef  int64_t ReadId;

typedef std::array<PosInRead, UPPER_KMER_FREQ> POSITIONS;
typedef std::array<ReadId,    UPPER_KMER_FREQ> READIDS;

typedef std::tuple<TKmer, int> KmerListEntry;
typedef std::vector<KmerListEntry> KmerList;

#define DEFAULT_THREAD_PER_TASK 4
#define MAX_THREAD_MEMORY_BOUNDED 8

/* 
 * Some explanations about these vars and recommended settings:
 * it's suggested to set number of MPI ranks per nodes as the count of NUMA nodes
 * under such setting, the memory bandwidth is not shared between MPI ranks
 * MAX_THREAD_MEMORY_BOUNDED_EXTREME is the count of memory controller per NUMA node ( for perlmutter, it is 2 )
 * For thr other vars, i don't have a good explanation for them currently
 */


struct KmerSeedStruct{
    TKmer kmer; 

    KmerSeedStruct(TKmer kmer) : kmer(kmer) {};
    KmerSeedStruct(const KmerSeedStruct& o) : kmer(o.kmer) {};
    // KmerSeedStruct(const KmerSeed& o) : kmer(std::get<0>(o)), readid(std::get<1>(o)), posinread(std::get<2>(o)) {};
    KmerSeedStruct(KmerSeedStruct&& o) :    // yfli: Not sure if it's appropriate to use std::move() here
        kmer(std::move(o.kmer)) {};
    KmerSeedStruct() {};

    int GetByte(int &i) const
    {
        return kmer.getByte(i);
    }

    bool operator<(const KmerSeedStruct& o) const
    {
        return kmer < o.kmer;
    }

    bool operator==(const KmerSeedStruct& o) const
    {
        return kmer == o.kmer;
    }

    bool operator!=(const KmerSeedStruct& o) const
    {
        return kmer != o.kmer;
    }

    KmerSeedStruct& operator=(const KmerSeedStruct& o)
    {
        kmer = o.kmer;
        return *this;
    }
};

typedef std::vector<std::vector<KmerSeedStruct>> KmerSeedBuckets;
typedef std::vector<std::vector<std::vector<KmerSeedStruct>>> KmerSeedVecs;


std::unique_ptr<KmerSeedBuckets> 
exchange_kmer(const DnaBuffer& myreads,
     int thr_per_task = DEFAULT_THREAD_PER_TASK,
     int max_thr_membounded = MAX_THREAD_MEMORY_BOUNDED);

std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBuckets>& recv_kmerseeds, 
     int thr_per_task = DEFAULT_THREAD_PER_TASK);

int GetKmerOwner(const TKmer& kmer, int ntasks);

void print_kmer_histogram(const KmerList& kmerlist);

struct KmerParserHandler
{
    int ntasks;
    std::vector<std::vector<KmerSeedStruct>>& kmerseeds;

    KmerParserHandler(std::vector<std::vector<KmerSeedStruct>>& kmerseeds) : ntasks(kmerseeds.size()), kmerseeds(kmerseeds) {}

    void operator()(const TKmer& kmer)
    {
        kmerseeds[GetKmerOwner(kmer, ntasks)].emplace_back(kmer);
    }
};


template <typename KmerHandler>
void ForeachKmerParallel(const DnaBuffer& myreads, std::vector<KmerHandler>& handlers, int nthreads)
{
    assert(nthreads > 0);

    /* cosidering the NUMA effect, we may want the vecs to be NUMA-local*/
    #pragma omp parallel for num_threads(nthreads) 
    for (size_t i = 0; i < myreads.size(); ++i)
    {
        int tid = omp_get_thread_num();

        if (myreads[i].size() < KMER_SIZE)
            continue;

        std::vector<TKmer> repmers = TKmer::GetRepKmers(myreads[i]);

        size_t j = 0;

        for (auto meritr = repmers.begin(); meritr != repmers.end(); ++meritr, ++j)
        {
            handlers[tid](*meritr);
        }
    }
}

#endif // KMEROPS_HPP