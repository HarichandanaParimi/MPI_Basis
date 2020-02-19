#include <mpi.h>
#include <unistd.h>
#include <iostream>

int main(int argc, char*argv[]) {
int rank , size;
char hostname[256];
   hostname[255] = '\0';
   gethostname(hostname, 1023);
MPI_Init(&argc,& argv);
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_Comm_size(MPI_COMM_WORLD, &size);
std::cout<<" i am process%"<<rank<<"out of."<<size<<"i am running on %s.";
MPI_Finalize();
}
