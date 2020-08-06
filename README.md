# distributed-branch-and-bound

Branch and Bound is an algorithm for Mixed Integer Linear Programming
Here, we have a developed a version of Branch and Bound which runs on parallel architectures and hence is sped up considerably.

The code for the Serial and Parallel versions is contained in two files respectively - Serial.cc and Parallel.cc

First of all, rename the Makefile_Serial or Makefile_Parallel to Makefile, depending upon the version you wish to compile.

Then

For Serial -

make build SOURCE=Serial.cc

./bin/Serial

For Parallel -

make build SOURCE=Parallel.cc

mpirun -np 8 ./bin/Parallel
