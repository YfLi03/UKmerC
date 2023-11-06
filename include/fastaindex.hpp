#ifndef FASTA_INDEX_H_
#define FASTA_INDEX_H_

#include "dnabuffer.hpp"

class FastaIndex
{
public:
    typedef struct { size_t len, pos, bases; } Record;

    FastaIndex(const std::string& fasta_fname);

    std::string get_fasta_fname() const { return fasta_fname; }
    std::string get_faidx_fname() const { return fasta_fname + ".fai"; }

    size_t getmyreadcount() const { return readcounts; }
    std::vector<size_t> getmyreadlens() const;

    const std::vector<Record>& getmyrecords() const { return myrecords; }

    DnaBuffer getmydna() const;
    void log(const DnaBuffer& buffer) const;

    static Record get_faidx_record(const std::string& line, std::string& name);


private:
    std::vector<Record> myrecords; /* records for the reads local processor is responsible for */
    std::vector<std::string> rootnames;
    uint64_t readcounts; /* number of reads assigned to each processor. Each processor gets a copy. |readcounts| == nprocs */
    std::string fasta_fname; /* FASTA file name */
};

#endif