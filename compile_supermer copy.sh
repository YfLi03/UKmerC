mkdir bin
make clean
make K=17 M=10 L=2 U=50 LOG=2 BATCH=100000 -j8
mv ukmerc bin/smpi_k17_m10_l2_u50
make clean
make K=31 M=17 L=2 U=50 LOG=2 BATCH=100000 -j8
mv ukmerc bin/mpi_k31_m17_l2_u50
make clean
make K=55 M=23 L=2 U=50 LOG=2 BATCH=100000 -j8
mv ukmerc bin/mpi_k55_m23_l2_u50
make clean
