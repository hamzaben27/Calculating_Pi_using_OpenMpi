// mpicc pi_omp_cal.c -lm -o piomp
// to run mpirun piomp
#include <stdio.h>
#include <time.h>
#include<mpi.h>
#include<math.h>
#include<omp.h>
#define f(A) (4.0/(1.0+A*A))
#define PI25DT 3.141592653589793238462643
const int n = 1 << 2 ; /* 1024 * 128 * 128 */

  const int WORKTAG=1;
  const int DIETAG=2 ;
int main ()
{
  int j ;
int fin = 40,debut = 14 ;
double width =  (fin-debut)/ n ;  
double  x1 = debut + j * width ;
double x2 = debut + (j+1)* width ;
int Ci = ((x1+x2)*width)/2 ; 
  int myrank ; 
  float fx = 4*(atan(fin)-atan(debut)); 
  MPI_Status status; 
  float PI_FIn = 4*fx;

  
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 

  //cacule ma valeur approx du pi
   
    if (myrank == 0) {

    double startwtime, endwtime;
    startwtime = MPI_Wtime();
    for(int i=1;i<n;i++)
    {
    MPI_Recv(&fx,1,MPI_FLOAT, 0, 10, MPI_COMM_WORLD, &status);
     
      printf(" pi %f ",fx);
      PI_FIn = PI_FIn+fx ;
    }
    printf("PI %f",PI_FIn);
    endwtime = MPI_Wtime();
    printf ( "Wall clock time = %f \n", endwtime-startwtime );
    }

     else {

      MPI_Send(&fx,1,MPI_FLOAT,myrank,WORKTAG,MPI_COMM_WORLD);
     }
    /* shut down mpi*/
    MPI_Finalize();

    return 0;
}

