#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "inverse.h"
#include "input.h"

void output(int l, int n, int m, double *A);
double norm(int n, int first_row, int last_row, double *B_loc, double *X_loc);
int min(int x, int y);
//int main1(int argc, char *argv[]);


int main(int argc, char *argv[]){ //argc - количество параметров командной строки, argv - массив строк - сами параметры, argv[0] всегда является именем программы
	
	double t;             // астрономическое время работы всего процесса
	
	double epsilon, epsilon_loc; // B - копия А, чтобы в конце посчитать норму невязки
	int i, n, m, k, p, curr, error=0, j=0;		// p - общее число процессов, curr - номер текущего процесса; порядок задания аргументов в командной строке такой же, но p - уже не параметр командной строки
	int first_row, last_row, max_rows, rows;
	double *A_loc, *X_loc, *B_loc, *A, *X; // чтобы посчитать норму невязки все равно понадобится весь массив
	
	MPI_Init(&argc, &argv);
	
	if ( (argc != 5) && (argc != 4) )  { //потому что имя программы тоже считается
		printf("Error! Incorrecr number of parameters\n");
		error = 1;
	}
	
	MPI_Comm_rank(MPI_COMM_WORLD, &curr);
	MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	n = atoi(argv[1]);
	m = atoi(argv[2]);
	k = atoi(argv[3]);
	
	
	
	if ( ( (k == 0) && (argc != 5) ) || ( (k != 0) && (argc != 4) ) || (p <= 0) ) {
		error = 1;	
	} 
	
	if (error == 1) {
		if (curr == 0) printf("Error! Incorrect parameters!\n");
		MPI_Finalize();
		return -1;
	}
	//printf("0\n");
	if (p > n)  p = n;
	
	first_row = n * curr;
	first_row /= p;
	
	last_row = n * (curr + 1);
	last_row = last_row / p - 1;
	
	rows = last_row - first_row + 1;
	max_rows = n/p + n%p;   
	
	A_loc = (double*)malloc(max_rows*n*sizeof(double)); // т.к будем пересылать между процессами, выделяем максимум 
	if (A_loc==NULL){
		error = 1;
		MPI_Bcast(&error, 1, MPI_INT, curr, MPI_COMM_WORLD);
	}
	//printf("1\n");
	if (error == 1) {
		if (curr == 0) printf("Not enough memory! \n");
		MPI_Finalize();
		return -1;
	}
	
	B_loc = (double*)malloc(max_rows*n*sizeof(double)); // копия А_loc
	if (B_loc==NULL){
		error = 1;
		MPI_Bcast(&error, 1, MPI_INT, curr, MPI_COMM_WORLD);
	}
	
	if (error == 1) {
		free(A_loc);
		if (curr == 0) printf("Not enough memory! \n");
		MPI_Finalize();
		return -1;
	}
	//printf("2\n");
	X_loc = (double*)malloc(max_rows*n*sizeof(double));
	if (X_loc==NULL){ 
		error = 1;
		MPI_Bcast(&error, 1, MPI_INT, curr, MPI_COMM_WORLD); 
	}
	
	if (error == 1) {
		free(A_loc);
		free(B_loc);
		if (curr == 0) printf("Not enough memory! \n");
		MPI_Finalize();
		return -1;
	}
	
	X = (double*)malloc(n*n*sizeof(double));
	if (X==NULL){ 
		error = 1;
		MPI_Bcast(&error, 1, MPI_INT, curr, MPI_COMM_WORLD); 
	}
	if (error == 1) {
		free(A_loc);
		free(B_loc);
		free(X_loc);
		if (curr == 0) printf("Not enough memory! \n");
		MPI_Finalize();
		return -1;
	}
	
	//printf("3\n");	
	if (k!=0) {   // заполняем A_loc
		matr_input(k, n, A_loc, first_row, last_row);
	}
	else if (curr == 0) {
		A = (double*)malloc(n*n*sizeof(double));
		if (A == NULL) {
			printf("Error! Not enough memory! \n");
			error = 1;
		}
		if (error == 0) {
			i = matr_file(n, argv[4], A);
			if (i == -1) {
				printf("The file cannot be opened! \n");
				free(A);
				error = 1;
			}
			if (i == -2) {
				printf("The file doesn't meet the requirements! \n");
				free(A);
				error = 1;
			}
		}
	} 
	
	if (error == 1) {
		free(A_loc);
		free(B_loc);
		free(X_loc);
		free(X);
		MPI_Finalize();
		return -1;
	}
	//printf("4\n");
	if ((k==0) && (curr == 0)) {
		MPI_Scatter(&A, n*n, MPI_DOUBLE, &A_loc, n*max_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD); // стр 276
	}
		
	for (i=0; i<n*rows; i++) {      // заполняем B_loc
		B_loc[i] = A_loc[i];
	}
	//printf("5\n");
	j = 0;	
	for (i=first_row*n; i<(last_row+1)*n; i++, j++) {     // заполняем X_loc
		//printf("%d\n", j);
		if ((i/n) == (i%n)) X_loc[j] = 1;
		else X_loc[j] = 0;
	}
	
	if (curr == 0) {
		printf("The matrix A : \n");
		output(rows, n, m, A_loc);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);  // как синхронайз
	t = MPI_Wtime();
	printf("6\n");
	if ( (inverse_matrix(n, A_loc, X_loc, curr, first_row, last_row, max_rows, p) ) != 0 ) {
		error = 1;
	}	
	printf("7\n");
	MPI_Barrier(MPI_COMM_WORLD); 
	t = MPI_Wtime() - t;
	
	if (error == 1) MPI_Bcast(&error, 1, MPI_INT, curr, MPI_COMM_WORLD);
	if (error == 1) {
		free(A_loc);
		free(B_loc);
		free(X_loc);
		if (curr == 0) printf("Error! The matrix is degenerate\n");
		MPI_Finalize();
		return -1;
	}
		
	
	
	if (curr == 0) {
		printf("The matrix A^(-1) : \n");
		output(rows, n, m, X_loc);
	}
	if (curr == 0) printf("Program execution time :  %lf \n", t);
	
	/*MPI_Allgather(&X_loc, n*max_rows, MPI_DOUBLE, &X, n*n, MPI_DOUBLE, MPI_COMM_WORLD);
	*epsilon_loc = norm(n, first_row, last_row, B_loc, X);
	MPI_Reduce(&epsilon_loc, &epsilon, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);  
	if (curr == 0) {
		*epsilon = sqrt(*epsilon);
		printf("The norm of the residual equals  %10.3e \n", *epsilon);
	}
	*/
	free(A_loc);
	free(B_loc);
	free(X_loc);
	if (k==0) free(A);

	MPI_Finalize();
	
	return 0;		
}


void output(int l, int n, int m, double *A){

	int i, counter_i = 0;
	int newl = min(l, m); // l - количество строк
	int newn = min(n, m); 
		
	for (i=0; i<newl*newn; i++){ //чтобы не выходил за рамки матрицы
		printf(" %10.3e", A[i+(n-newn)*counter_i]);
		if ( (i+1)%newn == 0 ) { //то есть указатель пришел на конец одной из строк выводимого среза матрицы
			counter_i += 1;
			printf("\n");
		}
	}
}


int min(int x, int y){
	if (x>=y) return y;
	else return x;
}

double norm(int n, int first_row, int last_row, double *B_loc, double *X) {
	int i, j, i_loc, k;
	double p, norm_loc = 0;
	
	for (i=first_row; i<last_row+1; i++) {
		i_loc = i-first_row;
		
		for (j=0; j<n; j++) {
			p = 0;
			for (k=0; k<n; k++) {
				p += B_loc[i_loc*n+j] * X[j*n+k];
			}
			if  (i == j) norm_loc += (p-1)*(p-1);
			else norm_loc += p*p;
		}
	}

	return norm_loc;
}


