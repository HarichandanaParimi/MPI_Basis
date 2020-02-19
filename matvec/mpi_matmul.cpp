#include <iostream>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <mpi.h>
#include <math.h>

using namespace std;

float genA (int row, int coloumn) {
  if (row > coloumn)
    return 1.0;
  else
    return 0.0;
}

float genx0 (int i) {
  return 1.0;
}

void checkx (int iter, long i, float xval) {
  if (iter == 1) {
    float shouldbe = i;
    if (fabs(xval/shouldbe) > 1.01 || fabs(xval/shouldbe) < .99 )
      cout<<"incorrect : x["<<i<<"] at iteration "<<iter<<" should be "<<shouldbe<<" not "<<xval<<endl;
  }

  if (iter == 2) {
    float shouldbe =(i-1)*i/2;
    if (fabs(xval/shouldbe) > 1.01 || fabs(xval/shouldbe) < .99)
      cout<<"incorrect : x["<<i<<"] at iteration "<<iter<<" should be "<<shouldbe<<" not "<<xval<<endl;
  }
}

//perform dense y=Ax on an n \times n matrix
void matmul(float*A, float*x, float*y, long n) {

  for (long row = 0; row<n; ++row) {
    float sum = 0;

    for (long coloumn = 0; coloumn<n ; ++coloumn) {
       sum += x[coloumn] * A[row*n+coloumn];
    }

       y[row] = sum;
  }
}

int main (int argc, char* argv[]) {

  if (argc < 3) {
    std::cout<<"usage: "<<argv[0]<<" <n> <iteration>"<<std::endl;
  }


  //initialize data
   bool check = true;

   long n = atol(argv[1]);

   int iter = atoi(argv[2]);

   MPI_Init(&argc,&argv);

   int w_rank,npsize;

   MPI_Comm_rank(MPI_COMM_WORLD,&w_rank);
   MPI_Comm_size(MPI_COMM_WORLD,&npsize);

   int p = sqrt(npsize);
   long e_div = n/p;

   int r_d = w_rank/p,c_d = w_rank%p;

   MPI_Comm rowwise_comm;
   int r_rank;
   MPI_Comm_split(MPI_COMM_WORLD, r_d, w_rank, &rowwise_comm);
   MPI_Comm_rank(rowwise_comm,&r_rank);

   MPI_Comm colwise_comm;
   int colrank;
   MPI_Comm_split(MPI_COMM_WORLD, c_d, w_rank, &colwise_comm);
   MPI_Comm_rank(colwise_comm,&colrank);

   float* A = new float[e_div*e_div];
   long rowwise_start = (r_d*e_div),colclock_start = (c_d*e_div);
   long r_end = rowwise_start+e_div,colend = colclock_start+e_div;

  for (long row = rowwise_start,rowwiseset=0; row<r_end; row++,rowwiseset++) {
    for (long coloumn= colclock_start,colwiseset=0; coloumn<colend; coloumn++,colwiseset++) {
      A[(rowwiseset*e_div)+colwiseset] = genA(row, coloumn);
    }
  }

  float* x = new float[e_div];

  for (long i=0; i<e_div; i++)
    {
     x[i] = genx0(i);
    }
  float* y = new float[e_div];
  for (long i=0; i<e_div; i++)
    y[i] = 0.0;

   double clock_start = MPI_Wtime();
   for (int k = 0; k<iter; k++)
   {
      matmul(A,x,y,e_div);
      MPI_Reduce(y,x,e_div,MPI_FLOAT,MPI_SUM,r_d,rowwise_comm);
      MPI_Bcast(x,e_div,MPI_FLOAT,c_d, colwise_comm);
	if (check){
	  for (long i = colclock_start,p=0; i<colend; ++i,++p){
	      checkx (k+1, i, x[p]);
             }
          }
  }

  if(w_rank == 0){
        double clock_end = MPI_Wtime();
        cerr<<clock_end-clock_start<<endl;
   }

  MPI_Comm_free(&rowwise_comm);
  MPI_Comm_free(&colwise_comm);

  delete[] A;
  delete[] x;
  delete[] y;

  MPI_Finalize();

  return 0;
}
