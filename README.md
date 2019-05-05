# Distributed-Branch-and-Bound

README file containing instructions to compile and run the code must be uploaded. 

The code for the Serial and Parallel versions is contained in two files respectively - Serial.cc and Parallel.cc

The code to compile and run serial file from the root directory.
First of all, rename the Makefile_Serial or Makefile_Parallel to Makefile, depending upon the version you wish to compile.

After that,
For Serial -
make build SOURCE=Serial.cc
./bin/Serial

For Parallel-
make build SOURCE=Parallel.cc
mpirun -np 8 ./bin/Parallel
