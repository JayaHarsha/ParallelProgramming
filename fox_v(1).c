
#include <mpi.h>
#define MAX 65536
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {int count;
float input[MAX];
#define siz(A) ((A)->count)
#define Entry(A,i,j) (*(((A)->input) + ((A)->count)*(i) + (j)))
} Mat_name;

typedef struct {int i ,j, row,col,rank;         
MPI_Comm  com, rcom,ccom;      
} mat_func;

void   mult_mat(Mat_name* mat1, Mat_name* mat2, Mat_name* mat3);
void   Mat_Generator(Mat_name* mat1, mat_func* matrix, int n);
Mat_name*  Mat_Sz_dec(int count);
void   create_mat(Mat_name* mat1);
MPI_Datatype     type_mat;
void   Return_output(Mat_name* mat1, mat_func* matrix, int n);
Mat_name*  temp_mat;
void   replicate(Mat_name* mat1);

int main(int argc, char* argv[]) {
	int rank, n, count;
	int array[10] = {100,200,300,400,500,600,700,800,900,1000};
	double entime,tottime,statime;
	mat_func matrix;
	Mat_name* mat1;
	Mat_name* mat2;
	Mat_name* mat3;
	void Setup_matrix(mat_func*  matrix);
	void multiplication(int n, mat_func* matrix, Mat_name* mat1,Mat_name* mat2, Mat_name* mat3);
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	Setup_matrix(&matrix);
	for(int loop = 0;loop < 10;loop++) {
	if (rank == 0) {
		n = array[loop];
		if(n % matrix.j != 0)
				exit(0);
	}
	statime = MPI_Wtime();
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	count = n/matrix.j;
	mat1 = Mat_Sz_dec(count);
	siz(mat1) = count;
	Mat_Generator(mat1, &matrix, n);
	/*Return_output(mat1, &matrix, n);*/
	mat2 = Mat_Sz_dec(count);
	siz(mat2) = count;
	Mat_Generator(mat2, &matrix, n);
	/*Return_output(mat2, &matrix, n);*/
	create_mat(mat1);
	temp_mat = Mat_Sz_dec(count);
	mat3 = Mat_Sz_dec(count);
	siz(mat3) = count;
	multiplication(n, &matrix, mat1, mat2, mat3);
	tottime = MPI_Wtime() - statime;
	/*Return_output(mat3, &matrix, n);*/
	printf("\n");
	if(rank == 0) printf("total time taken for the %d processors for %d is %lf\n",matrix.i,array[loop],tottime);
	}
	MPI_Finalize();
	
	return 0;
}  

void Setup_matrix(mat_func*  matrix ) {
	int my_rank, dims[2],periods[2],ords[2], remain_dims[2];
	MPI_Comm_size(MPI_COMM_WORLD, &(matrix->i));
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	matrix->j = (int) sqrt((double) matrix->i);
	dims[0] = dims[1] = matrix->j;
	periods[0] = periods[1] = 1;
	MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &(matrix->com));
	MPI_Comm_rank(matrix->com, &(matrix->rank));
	MPI_Cart_coords(matrix->com, matrix->rank, 2, ords);
	matrix->row = ords[0];
	matrix->col = ords[1];
	remain_dims[0] = 0; 
	remain_dims[1] = 1;
	MPI_Cart_sub(matrix->com, remain_dims, &(matrix->rcom));
	remain_dims[0] = 1; 
	remain_dims[1] = 0;
	MPI_Cart_sub(matrix->com, remain_dims, &(matrix->ccom));
} 

void multiplication(int n, mat_func* matrix, Mat_name*  mat1  ,Mat_name*  mat2  ,Mat_name*  mat3  ) {
	Mat_name*  new_m; 
	int i, src, count, s, d;
	MPI_Status status;
	count = n/matrix->j;
	replicate(mat3);
	s = (matrix->row + 1) % matrix->j;
	d = (matrix->row + matrix->j - 1) % matrix->j;  
	new_m = Mat_Sz_dec(count);
	for (i = 0; i < matrix->j; i++) {
		src = (matrix->row + i) % matrix->j;
		if (src == matrix->col) {
			MPI_Bcast(mat1, 1, type_mat,src, matrix->rcom);
			mult_mat(mat1, mat2, mat3);
		} else {
			MPI_Bcast(new_m, 1, type_mat,src, matrix->rcom);
			mult_mat(new_m, mat2, mat3);
		}
		MPI_Sendrecv_replace(mat2, 1, type_mat,d, 0, s, 0, matrix->ccom, &status);
	} 
} 

Mat_name* Mat_Sz_dec(int local_order) {
	Mat_name* temp;  
	temp = (Mat_name*) malloc(sizeof(Mat_name));
	return temp;
} 

void Mat_Generator(Mat_name*  mat1, mat_func* matrix,int n) {
	int mr, mc, matr, matc, d, ords[2];
	float* temp;
	MPI_Status status;	
	if (matrix->rank == 0) {
		temp = (float*) malloc(siz(mat1)*sizeof(float));
		fflush(stdout);
		for(mr = 0;  mr < n; mr++) {
			matr = mr/siz(mat1);
			ords[0] = matr;
			for (matc = 0; matc < matrix->j; matc++) {
				ords[1] = matc;
				MPI_Cart_rank(matrix->com, ords, &d);
				if (d == 0) {
					for(mc = 0; mc < siz(mat1); mc++)
						(*((mat1->input)+mr*siz(mat1)+mc)) = (((rand()/((float)RAND_MAX)) - 0.5)*2);
				} else {
					for(mc = 0; mc < siz(mat1); mc++)
						(*(temp + mc)) = (((rand()/((float)RAND_MAX)) - 0.5)*2);
					MPI_Send(temp, siz(mat1), MPI_FLOAT, d, 0,matrix->com);
				}
			}
		}
		free(temp);
	} else {
		for (mr = 0; mr < siz(mat1); mr++) 
			MPI_Recv(&Entry(mat1, mr, 0), siz(mat1), MPI_FLOAT, 0, 0, matrix->com, &status);
	}
					 
} 

/*void Return_output( Mat_name*  mat1,mat_func* matrix,int n) {
	int mr, mc, matr, matc,s, ords[2];
	float* temp;
	MPI_Status status;
	if (matrix->rank == 0) {
		temp = (float*) malloc(siz(mat1)*sizeof(float));
		for (mr = 0;  mr < n; mr++) {
			matr = mr/siz(mat1);
			ords[0] = matr;
			for (matc = 0; matc < matrix->j; matc++) {
				ords[1] = matc;
				MPI_Cart_rank(matrix->com, ords, &s);
				if (s == 0) {
					for(mc = 0; mc < siz(mat1); mc++)
						printf("%4.1f ", Entry(mat1, mr, mc));
				} else {
					MPI_Recv(temp, siz(mat1), MPI_FLOAT, s, 0,matrix->com, &status);
					for(mc = 0; mc < siz(mat1); mc++)
						printf("%4.1f ", temp[mc]);
				}
			}
			printf("\n");
		}
		free(temp);
	} else {
		for (mr = 0; mr < siz(mat1); mr++) 
			MPI_Send(&Entry(mat1, mr, 0), siz(mat1), MPI_FLOAT, 0, 0, matrix->com);
	}
					 
}  */

void replicate(Mat_name*  mat1  ) {
	int i, j;
	for (i = 0; i < siz(mat1); i++)
		for (j = 0; j < siz(mat1); j++)
			Entry(mat1,i,j) = 0.0;
} 

void create_mat(Mat_name*  mat1 ) {
	MPI_Datatype  type_temp;
	int lens[2];
	MPI_Aint vecs[2], stad,address;
	MPI_Datatype  list[2];
	MPI_Type_contiguous(siz(mat1)*siz(mat1), MPI_FLOAT, &type_temp);
	lens[0] = lens[1] = 1;
	list[0] = MPI_INT;
	list[1] = type_temp;
	MPI_Get_address(mat1, &stad);
	MPI_Get_address(&(mat1->count), &address);
	vecs[0] = address - stad;
	MPI_Get_address(mat1->input, &address);
	vecs[1] = address - stad;
	MPI_Type_create_struct(2, lens, vecs,list, &type_mat);
	MPI_Type_commit(&type_mat); 
} 

void mult_mat(Mat_name*  mat1,Mat_name*  mat2, Mat_name*  mat3) {
	int a, b, c;
	for (a = 0; a < siz(mat1); a++)
		for (b = 0; b < siz(mat1); b++)
			for (c = 0; c < siz(mat2); c++)
				Entry(mat3,a,b) = Entry(mat3,a,b) + Entry(mat1,a,c)*Entry(mat2,c,b);
}  