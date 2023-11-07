#include "kmerops.hpp"
#include "logger.hpp"
#include "dnaseq.hpp"
#include "timer.hpp"
#include <cstring>
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <cmath>

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

std::unique_ptr<KmerSeedBucket>
extract_kmer(const DnaBuffer& myreads)
{
    Logger logger;
    std::ostringstream rootlog;

    size_t numreads = myreads.size();

    KmerSeedBucket* kmerseeds = new KmerSeedBucket;
    KmerParserHandler handler(*kmerseeds);

    ForeachKmer(myreads, handler);

    return std::unique_ptr<KmerSeedBucket>(kmerseeds);
}


std::unique_ptr<KmerList>
filter_kmer(std::unique_ptr<KmerSeedBucket>& kmerseeds)
{
    Timer timer;
    Logger logger;
    std::ostringstream rootlog;

    std::sort(kmerseeds->begin(), kmerseeds->end());

    uint64_t valid_kmer = 0;
    KmerList* kmerlist = new KmerList();
    
    TKmer last_mer = (*kmerseeds)[0].kmer;
    uint64_t cur_kmer_cnt = 1;

    for(size_t idx = 1; idx < (*kmerseeds).size(); idx++) {
        TKmer cur_mer = (*kmerseeds)[idx].kmer;
        if (cur_mer == last_mer) {
            cur_kmer_cnt++;
        } else {
            if (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ) {

                kmerlist->push_back(KmerListEntry());
                KmerListEntry& entry    = kmerlist->back();
                TKmer& kmer             = std::get<0>(entry);
                int& count              = std::get<1>(entry);

                count = cur_kmer_cnt;
                kmer = last_mer;
                valid_kmer++;
            }

            cur_kmer_cnt = 1;
            last_mer = cur_mer;
        }
    }

    /* deal with the last kmer */
    if (cur_kmer_cnt >= LOWER_KMER_FREQ && cur_kmer_cnt <= UPPER_KMER_FREQ) {
        kmerlist->push_back(KmerListEntry());
        KmerListEntry& entry         = kmerlist->back();
        TKmer& kmer             = std::get<0>(entry);
        int& count              = std::get<1>(entry);

        count = cur_kmer_cnt;
        kmer = last_mer;
        valid_kmer++;
    }
    
    #if LOG_LEVEL >= 2
    logger() << valid_kmer << " ";
    logger.flush("Valid kmer total:");
    #endif

    return std::unique_ptr<KmerList>(kmerlist);
}