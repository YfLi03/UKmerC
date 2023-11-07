#include<iostream>
#include <iomanip>
#include <sstream>

#include "logger.hpp"
#include "timer.hpp"
#include "fastaindex.hpp"
#include "dnabuffer.hpp"
#include "dnaseq.hpp"
#include "kmerops.hpp"

std::string fasta_fname;

int main(int argc, char **argv){

    Logger log;
    Timer timer;
    std::ostringstream ss;

    if (argc < 2){
        std::cerr << "Usage: " << argv[0] << " <fasta file>" << std::endl;
        exit(1);
    }

    fasta_fname = argv[1];

    /* read index */
    timer.start();
    FastaIndex index(fasta_fname);
    ss << "reading " << std::quoted(index.get_faidx_fname()) << " and scattering to all MPI tasks";
    timer.stop_and_log(ss.str().c_str());
    ss.clear(); ss.str("");

    /* read fasta and encode */
    timer.start();
    DnaBuffer mydna = index.getmydna();
    ss << "reading and 2-bit encoding " << std::quoted(index.get_fasta_fname()) << " sequences in parallel";
    timer.stop_and_log(ss.str().c_str());
    ss.clear(); ss.str("");

    /* extract kmers from sequences*/
    timer.start();
    auto bucket = extract_kmer(mydna);
    timer.stop_and_log("extract_kmer");

    /* filter valid kmers */
    timer.start();
    auto kmerlist = filter_kmer(bucket);
    timer.stop_and_log("filter_kmer");

    /* print result */
    print_kmer_histogram(*kmerlist);

    return 0;
}