#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(void){

        int a,b,c,i,count;
        int n = 32;
        int comm_sz;
        int my_rank;
        double statime,entime;
        double *A;
        double *B;
        double *product;
        double p;
        int rowindex,colindex;
        int append=0;
        MPI_Status status;
        MPI_Init(NULL,NULL);
        MPI_Comm_size(MPI_COMM_WORLD,&comm_sz);
        MPI_Comm_rank(MPI_COMM_WORLD,&my_rank);

        A = malloc(n*n*sizeof(double));
        B = malloc(n*n*sizeof(double));
        product = malloc(n*n*sizeof(double));

        srand(0);
        if(my_rank == 0){
                for(i = 0; i<n*n;i++){
                        A[i] = 0 + (1 - 0) * (rand() /((double) RAND_MAX));

                        B[i] = 0 + (1 - 0) * (rand() /((double) RAND_MAX));
                }



        }
        count = n*n/comm_sz;
        double *temprow = malloc(count*sizeof(double));
        double *tempcol = malloc(count*sizeof(double));
        statime = MPI_Wtime();
        MPI_Scatter(A,count,MPI_DOUBLE,temprow,count,MPI_DOUBLE,0,MPI_COMM_WORLD);

        MPI_Scatter(B,count,MPI_DOUBLE,tempcol,count,MPI_DOUBLE,0,MPI_COMM_WORLD);

        double *temprow1 = malloc(n*sizeof(double));
        double *tempcol1 = malloc(n*sizeof(double));
        for(a = 2*my_rank;a < (n/comm_sz)+2*my_rank;a++){
		 for(b = 2*my_rank;b < (n/comm_sz) +2*my_rank ;b++){

                        temprow1[b] = temprow[b +n*a];

                        tempcol1[b] = tempcol[b +n*a];
                        p = 0;
                        for(c = 0;c < n;c++){
                                p = p + temprow1[c]*tempcol1[c];

                        }
                product[append] = p;
                append++;
                }
        }

		for(i =1;i<comm_sz;i++){
        if(my_rank == 0){
        MPI_Send(temprow,count,MPI_DOUBLE,my_rank+1,0,MPI_COMM_WORLD);
        MPI_Recv(temprow,count,MPI_DOUBLE,comm_sz-1,0,MPI_COMM_WORLD,&status);
        }
        else if(my_rank == comm_sz-1 ){
         MPI_Send(temprow,count,MPI_DOUBLE,0,0,MPI_COMM_WORLD);
         MPI_Recv(temprow,count,MPI_DOUBLE,my_rank -1,0,MPI_COMM_WORLD,&status);

        }
        else{
         MPI_Send(temprow,count,MPI_DOUBLE,my_rank+1,0,MPI_COMM_WORLD);
         MPI_Recv(temprow,count,MPI_DOUBLE,my_rank-1,0,MPI_COMM_WORLD,&status);

        }

         for(a = 0;a < n/comm_sz;a++){
                 for(b = 0;b < n/comm_sz;b++){

                        temprow1[b] = temprow[b +n*a];

                        tempcol1[b] = tempcol[b +n*a];
                        p = 0;
                        for(c = 0;c < n;c++){
                                p = p + temprow1[c]*tempcol1[c];

                        }
                product[append] = p;
                append++;
                }
        }

}

  double* final_product = malloc(n*n*comm_sz*sizeof(double));
        double* matrixproduct = malloc(n*n*sizeof(double));
entime = MPI_Wtime() - statime;
        if(my_rank == 0)
        {
        printf("Time for %d processors is %lf\n",comm_sz, entime);
        }
        MPI_Finalize();
        return 0;

}
