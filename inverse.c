#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "inverse.h"

void output(int l, int n, int m, double *A);
double mod(double c);

int inverse_matrix(int n, double *A_loc, double *X_loc, int curr, int first_row, int last_row, int max_rows, int p){ 
	double a, tmp, a_loc, b;					    
	int i, i_loc, j, k, m, m_loc, proc, flag, k_loc;  
	double n_A, norm_loc;
	double *row_to_send;
	MPI_Status st;
	
	row_to_send = (double*)malloc(2*n*sizeof(double));
	
	norm_loc = norma(n, A_loc, last_row-first_row+1); // посчитали норму всей матрицы и передали всем процессам
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allreduce(&norm_loc, &n_A, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	printf("0\n");
	
	flag = 0;
	
	for (i = 0; i<n; i++){ // (n-i) - ранг подматрицы
		
		if (i <= last_row) {	     // упорядочиваем матрицу
			if ( first_row <= i )  {
				m = i; // номер строки, содержащей максимальный элемент столбца i в подматрице
			}
			else m = first_row;      // но в большой матрице А
			
			k_loc = m - first_row;   // а это уже для A_loc
			a_loc = A_loc[k_loc * n + i];			      // максимальный элемент в столбце i
			printf("i = %d,  k_loc = %d\n",i, k_loc);
			for (j=k_loc+1; j<last_row-first_row+1; j++){                      // находим максимальный элемент в столбце подматрицы
				printf("a_loc = %lf\n", a_loc);
				
				if ( fabs(A_loc[j * n + i]) > a_loc )  {
					k_loc = j; 
					printf("i = %d,  k_loc = %d\n",i, k_loc);
					a_loc = A_loc[j * n + i];
				}		
			}
			
		}
		MPI_Barrier(MPI_COMM_WORLD);  
		MPI_Allreduce(&a_loc, &a, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
		// упорядочили
		printf("%lf\n %lf\n", a_loc, a);
		
		
		MPI_Barrier(MPI_COMM_WORLD);  
		if (a_loc >= a) {
			proc = curr;
			k = k_loc;
			MPI_Bcast(&k, 1, MPI_INT, curr, MPI_COMM_WORLD);
			MPI_Bcast(&proc, 1, MPI_INT, curr, MPI_COMM_WORLD);
		}
		
		printf("k = %d\n", k);
		
		if (curr == proc) {
			//printf("a = %lf \n, n_A = %lf\n", a, *n_A);
			if (fabs(a)/n_A < 2.20e-15) {
				flag = 1;
				MPI_Bcast(&flag, 1, MPI_INT, curr, MPI_COMM_WORLD);
			}
		}
		if (flag == 1) return -1;
		//printf("1\n");
		
		if (i <= last_row) {
			if ((proc == curr) && (i >= first_row)) { //  то есть максимум лежит в этом процессе
				b = 1./a; 	 
				i_loc = i - first_row;
				for (j=0; j<n; j++){ // переставили
						tmp = A_loc[i_loc * n + j];
						A_loc[i_loc * n + j] = A_loc[k * n + j];
						A_loc[k * n + j] = tmp;
						tmp = X_loc[i_loc * n + j];
						X_loc[i_loc * n + j] = X_loc[k * n + j];
						X_loc[k * n + j] = tmp;
				}
				for (j=0; j<n; j++) { // сразу нормировали
					A_loc[i_loc * n + j] = A_loc[i_loc * n + j] * b;
					X_loc[i_loc * n + j] = X_loc[i_loc * n + j] * b;
				}
				printf("X = \n");
				output(last_row-first_row+1, n, 5, X_loc);
				printf("A = \n");
				output(last_row-first_row+1, n, 5, A_loc);
			}
			
			if ( (proc != curr) && (i >= first_row) ) { // то есть мы должны i-ую строчку поменять с k-той, которая из другого процесса
				i_loc = i - first_row ;
				for (j=0; j<n; j++) {
					row_to_send[j] = A_loc[i_loc*n + j];
					row_to_send[j+n] = X_loc[i_loc*n + j];
				}
			}
			
			if ( (proc == curr) && (first_row > i) ) {
				b = 1./a ; 
				for (j=0; j<n; j++) {
					row_to_send[j] = A_loc[k * n + j]*b;
					row_to_send[j+n] = X_loc[k * n + j]*b;	
				}
			}
		}
		printf("2\n");
		MPI_Barrier(MPI_COMM_WORLD);
		// СОМНЕВАЮСЬ В СЛЕДУЮЩЕЙ СТРОЧКЕ
		// Вроде как можно заблочить только 2 процесса , но дешевле ли это, чем то, что дальше?
		if ( (proc != curr) && (i <= last_row) && (i >= first_row) ) {
			MPI_Sendrecv_replace(&row_to_send, 2*n, MPI_DOUBLE, proc, 0, curr, 0, MPI_COMM_WORLD, &st);
		}	  	
		
		if ( (curr != proc) && (i <= last_row) && (i >= first_row) ) {
			i_loc = i - first_row ;
			for (j=0; j<n; j++) {
				A_loc[i_loc*n + j] = row_to_send[j] ;
				X_loc[i_loc*n + j] = row_to_send[j+n];
			}
		}
		
		if ( (curr == proc) && (i > last_row) ) {
			for (j=0; j<n; j++) {
				A_loc[k*n + j] = row_to_send[j] ;
				X_loc[k*n + j] = row_to_send[j+n];
			}
		}
		
		// Упорядочили
		
		MPI_Barrier(MPI_COMM_WORLD);
		if ((first_row <= i) && (i <= last_row)) {
			for (j=0; j<n; j++) {
				row_to_send[j] = A_loc[i_loc*n + j];
				row_to_send[j+n] = X_loc[i_loc*n + j];
			}
			MPI_Bcast(&row_to_send, 2*n, MPI_DOUBLE, curr, MPI_COMM_WORLD);
		}
		
		 
		for (m=first_row; m<last_row+1; m++){                  // вычитаем из всех остальных с подходящим коэффициентом
			if (m != i){            // именно из остальных
				m_loc = m-first_row;
				a = A_loc[m_loc * n + i];
				for (j=n-1; j>-1; j--){  // чтобы было что вычитать(иначе первый элемент обнулится сразу)
					A_loc[m_loc * n + j] -= a * row_to_send[j];  // не делим на iый элемент indi[i] строчки, потому что он равен 1 из предыдущего шага
					X_loc[m_loc * n + j] -= a * row_to_send[j+n];
				}
			}
		}	
		
		MPI_Barrier(MPI_COMM_WORLD);	

		printf("X = \n");
		output(n, n, 5, X_loc);
		printf("A = \n");
		output(n, n, 5, A_loc);
	}
	
	return 0;

}


double norma(int n, double *A_loc, int rows) {
	int i, j;
	double sum_row, norm = 0;
	
	for (i=0; i<rows; i++) {
		sum_row = 0;
		for (j=0; j<n; j++){
			sum_row += mod(A_loc[i*n+j]);
		}
		if (sum_row >= norm) norm = sum_row;
	} 
	
	return norm;
}

double mod(double c){
	if (c >= 0) return c;
	else return (-c);
}
