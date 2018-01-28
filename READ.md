This package contains the codes for solving the heat equation using stencil based approach.
This code is written using MPI's one sided interface  for shared-memory configuration.

The usage are 
                     
mpirun -np 4 ./main.exe 10 20 2 2

arrgc=4
argv[1] n=10 ;   Number of grids in one direction.
argv[2] iter=10; Number of iteration for updates.
argv[3] px=2; Number of processes in X-direction
argv[4] py=2; Number of processes in Y-direction
