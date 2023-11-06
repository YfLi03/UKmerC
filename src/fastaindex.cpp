#include "fastaindex.hpp"
#include "logger.hpp"
#include "timer.hpp"
#include <cstring>
#include <iterator>
#include <algorithm>
#include <functional>
#include <numeric>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <memory>
#include <cassert>

using Record = typename FastaIndex::Record;

Record FastaIndex::get_faidx_record(const std::string& line, std::string& name)
{
    /*
     * Read a line from a FASTA index file into a record object.
     */
    Record record;
    std::istringstream(line) >> name >> record.len >> record.pos >> record.bases;
    return record;
}

FastaIndex::FastaIndex(const std::string& fasta_fname) :  fasta_fname(fasta_fname)
{
    /*
     * Root processor responsible for reading and parsing FASTA
     * index file "{fasta_fname}.fai" into one record per sequence.
     */

    std::string line, name;
    std::ifstream filestream(get_faidx_fname());

    while (std::getline(filestream, line))
    {
        myrecords.push_back(get_faidx_record(line, name));
        rootnames.push_back(name);
    }

    filestream.close();

    readcounts = myrecords.size();

    #if LOG_LEVEL >= 2
    Logger logger;
    logger() << " Reading index for " << readcounts << " sequences. ";
    logger.flush("Fasta index construction:");
    #endif
}

std::vector<size_t> FastaIndex::getmyreadlens() const
{
    /*
     *  Because we store sequence information using "array of structs"
     *  instead of "structs of arrays", if we want a linear array
     *  of read lengths we have to unpack them.
     */
    std::vector<size_t> readlens(getmyreadcount());
    std::transform(myrecords.cbegin(), myrecords.cend(), readlens.begin(), [](const auto& record) { return record.len; });
    return readlens;
}

DnaBuffer FastaIndex::getmydna() const
{
    /*
     * Allocate local sequence buffer.
     */
    auto readlens = getmyreadlens(); /* vector of local read lengths */
    size_t bufsize = DnaBuffer::computebufsize(readlens); /* minimum number of bytes needed to 2-bit encode all the local reads */
    DnaBuffer dnabuf(bufsize); /* initialize dnabuf by allocating @bufsize bytes */
    size_t numreads = readlens.size(); /* number of local reads */

    uint64_t startpos; /* the FASTA position that starts my local chunk of reads */
    uint64_t endpos; /* the FASTA position that ends my local chunk of reads (exclusive) */
    uint64_t filesize; /* the total size of the FASTA */
    uint64_t readbufsize; /* endpos - startpos */

    std::ifstream file(get_fasta_fname().c_str(), std::ios::binary); /* FASTA file handle */

    filesize = file.seekg(0, std::ios::end).tellg(); /* get the size of the FASTA file */

    startpos = myrecords.front().pos;
    endpos = myrecords.back().pos + myrecords.back().len + (myrecords.back().len / myrecords.back().bases);
    if (endpos > filesize) endpos = filesize;

    readbufsize = endpos - startpos;
    std::unique_ptr<char[]> readbuf(new char[readbufsize]);

    file.seekg(startpos, std::ios::beg);
    file.read(&readbuf[0], readbufsize);
    file.close();

    size_t totbases = std::accumulate(readlens.begin(), readlens.end(), static_cast<size_t>(0), std::plus<size_t>{});
    size_t maxlen = *std::max_element(readlens.begin(), readlens.end());

    /*
     * ASCII sequences are first read into this temporary char buffer
     * before they are compressed into the sequence buffer.
     */
    std::unique_ptr<char[]> tmpbuf(new char[maxlen]);

    #if LOG_LEVEL >= 2
    double st = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    #endif

    /*
     * Go through each local FASTA record.
     */
    for (auto itr = myrecords.cbegin(); itr != myrecords.cend(); ++itr)
    {
        size_t locpos = 0;
        ptrdiff_t chunkpos = itr->pos - startpos;
        ptrdiff_t remain = itr->len;
        char *writeptr = &tmpbuf[0];

        /*
         * Read ASCII FASTA sequence into the temoprary buffer.
         */
        while (remain > 0)
        {
            size_t cnt = std::min(itr->bases, static_cast<size_t>(remain));
            std::memcpy(writeptr, &readbuf[chunkpos + locpos], cnt);
            writeptr += cnt;
            remain -= cnt;
            locpos += (cnt+1);
        }

        /*
         * DnaBuffer automatically 2-bit encodes the ASCII sequence
         * and pushes it onto its local stack of sequences.
         */
        dnabuf.push_back(&tmpbuf[0], itr->len);
    }


    #if LOG_LEVEL >= 2
    double ed = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    double elapsed = ed - st;
    double mbspersecond = (totbases / 1048576.0) / elapsed;
    Logger logger;
    logger() << std::fixed << std::setprecision(2) << mbspersecond << " Mbs/second";
    logger.flush("FASTA parsing rates (DnaBuffer):");
    #endif

    return dnabuf;
}
