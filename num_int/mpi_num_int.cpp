#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <mpi.h>

using namespace std;

#ifdef __cplusplus
extern "C" {
#endif

float f1(float x, int intensity);
float f2(float x, int intensity);
float f3(float x, int intensity);
float f4(float x, int intensity);

#ifdef __cplusplus
}
#endif


int main (int argc, char* argv[]) {

  if (argc < 6) {
    std::cerr<<"usage: "<<argv[0]<<" <functionid> <a> <b> <n> <intensity>"<<std::endl;
    return -1;
  }

  MPI_Init(&argc,&argv);
  int functionid = atoi(argv[1]);
  float a = atof(argv[2]);
  float b = atof(argv[3]);
  int n = atoi(argv[4]);
  int intensity = atoi(argv[5]);
  int rank,p;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&p);
  float (*function)(float,int);
  float temp = (b-a)/n;
  int loopstart,loopend,partition;
  double clock_start = MPI_Wtime();
  partition = n/p;
  loopstart = rank*partition; loopend = (rank+1)*partition;
  if(rank == p-1)
	  loopend = n;
  switch(functionid)
  {
     case 1 : function = &f1;

                break;

     case 2 : function = &f2;
                break;

     case 3 : function = &f3;
                break;

     case 4 : function = &f4;
                break;

     default :  cerr<<"Functionid does not exists, Enter 1, 2, 3 or 4" <<endl;
                return -1;

   };
   float Partial_NumInt= 0.0,NumInt = 0.0;
   for(int i = loopstart ; i<loopend ;i++)
        {
           float x = a + ( (i + 0.5) *temp );
	       Partial_NumInt += (*function)(x,intensity);
        }
   MPI_Reduce(&Partial_NumInt,&NumInt,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
   if(rank == 0)
   {
	   NumInt = NumInt * temp;
	   cout<<NumInt<<endl;
           double clock_end = MPI_Wtime();
           cerr<<clock_end-clock_start<<endl;
   }
   MPI_Finalize();

  return 0;
}
