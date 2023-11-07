#ifndef KMEROPS_HPP
#define KMEROPS_HPP

#include "kmer.hpp"
#include "timer.hpp"
#include "dnaseq.hpp"
#include "dnabuffer.hpp"
#include "logger.hpp"

typedef std::tuple<TKmer, int> KmerListEntry;
typedef std::vector<KmerListEntry> KmerList;

struct KmerSeedStruct{
    TKmer kmer; 

    KmerSeedStruct(TKmer kmer) : kmer(kmer) {};
    KmerSeedStruct(const KmerSeedStruct& o) : kmer(o.kmer) {};
    KmerSeedStruct(KmerSeedStruct&& o) : kmer(std::move(o.kmer)) {};
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

typedef std::vector<KmerSeedStruct> KmerSeedBucket;


std::unique_ptr<KmerSeedBucket> 
extract_kmer(const DnaBuffer& myreads);

std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBucket>& kmerseeds);

int GetKmerOwner(const TKmer& kmer, int ntasks);

void print_kmer_histogram(const KmerList& kmerlist);

struct KmerParserHandler
{
    std::vector<KmerSeedStruct>& kmerseeds;

    KmerParserHandler(std::vector<KmerSeedStruct>& kmerseeds) : kmerseeds(kmerseeds) {}

    void operator()(const TKmer& kmer)
    {
        kmerseeds.emplace_back(kmer);
    }
};


template <typename KmerHandler>
void ForeachKmer(const DnaBuffer& myreads, KmerHandler& handler)
{
    for (size_t i = 0; i < myreads.size(); ++i)
    {

        if (myreads[i].size() < KMER_SIZE)
            continue;

        std::vector<TKmer> repmers = TKmer::GetRepKmers(myreads[i]);

        for (auto meritr = repmers.begin(); meritr != repmers.end(); ++meritr)
        {
            handler(*meritr);
        }
    }
}

#endif // KMEROPS_HPP