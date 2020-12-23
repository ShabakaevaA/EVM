#include <stdio.h>
#include <stdlib.h>
#include "input.h"

void matr_input(int k, int n, double *A, int first_row, int last_row) {
	int i, j=0;
	
	for (i = first_row * n ; i < (last_row+1) * n; i++, j++){
		A[j] = f(k, n, i/n + 1, i%n + 1);  // a[i][j] = A[i*n+j] - как для элементов матрицы, то есть i>= 1, j>= 1;  int/int = int
	}
}


double f(int k, int n, int i, int j) 
{
	if (k == 1) return ( n-max(i, j)+1 );
	if (k == 2) return ( max(i, j) );
	if (k == 3) return ( max(i-j, j-i) );
	if (k == 4) return ( 1./(i+j-1));
	if (k == 5) {
		if (j == n) return (0.) ;
		else return ( n-max(i, j)+1 );
	}	
}


int matr_file(int n, char *filename, double *A){
	FILE *IN ;
	int i;
	
	IN = fopen(filename, "r");
	if (IN == NULL) {
		return -1;
	}
	
	for (i=0; i<n*n; i++) {
		if (fscanf(IN, "%lf", &A[i]) != 1) {
			fclose(IN);
			return -2;
		}	
	}
	
	fclose(IN);
	return 0;
}

int max(int x, int y) {
	if ( x >= y ) return x;
	else return y;
}
		
