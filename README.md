## ukmerc

A simple serial demo for K-mer counting is provided here. It is suggested to use this as a framework to develop your parallel code. Utilities such as timer and logger are provided. Some currently unused functions provided, such as hash functions or `GetKmerOwner`, may come in handy.

The `extract_kmer` and `filter_kmer` function in `kmerops.cpp` is taking the longest time currently. Try to parallelize these two functions. You don't actually need to modify other code files. You can also introduce additional code files if necessary, just remember to modify the Makefile if you do that. 


To compile, use `make` in the project directory. For example, you can use

```sh
make K=51 L=15 U=40 -j8
```

where K indicates the length of K-mer, L indicates the lower bound for K-mer counting, and U indicates the upper bound for K-mer counting. Other compile parameters include LOG level and DEBUG level.


To run on a cluster (such as NERSC's perlmutter), use

```sh
srun -C cpu -N 1 -n 1 -c 128 --cpu_bind=cores ./ukmerc PATH_TO_YOUR_DATASET
```

Further explanation to the runtime parameters can be found at NERSC's documentation. Our binary application is only taking in one parameter, which is the path to the fasta dataset.


Note that your fasta dataset should be indexed in advance. To do this on perlmutter, execute the following commands:

```sh
module load cpu
module load spack
spack install samtools
spack load samtools

samtools faidx  PATH_TO_YOUR_DATASET        # this is the index step

module load gpu                             # restore the default environment
```