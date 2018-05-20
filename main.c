/*
 * Author: James Wu
 * Date: 19 May 2018
 * Purpose: To implement a Gauss-Jordan Elimination algorithm
 *          for an arbitraily sized m x n matrix A
 */

#define _CRT_SECURE_NO_WARNINGS

/* Header files */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Constants */
#define FALSE 0
#define TRUE 1
#define m 4
#define n 4
#define EPSILON 1e-6

#define FILENAME "rref.tex"
#define TEXAUTHOR "GJE2"
#define TEXDATE "19 May 2018"

/* Function Prototypes */
void rref(double A[m][n]);
void swap(double A[m][n], int a, int b);
void mult(double A[m][n], int a, double c);
void add(double A[m][n], int a, int b, double c);
void preambles(FILE* myFile);
void endambles(FILE* myFile);
void texmatrix(FILE* myFile, double A[m][n]);

int main(void) {

	/* INIT m x n matrix */										
	double A[m][n] = { {2, -5, -3, 16},
	                   {5, -6, 6, -13},
					   {-2, -3, 6, 10},
					   {23, -19, -33, 27} };

	printf("Writing to %s...\n", FILENAME);

	rref(A);

	printf("Finished writing to %s!\n", FILENAME);

	system("PAUSE");
	return 0;

}

/*
 * The master Gauss-Jordan Elimination algorithm.
 * Parameters: A (double[m][n]) - the m x n matrix to be put into rref
 */
void rref(double A[m][n]) {

	FILE* myFile = fopen(FILENAME, "w");

	int anchor = 0; // row # that next anchor should be. Note: anchor = pivot
	int anchorExists; // true if anchor for given column exists, false 

	preambles(myFile);

	fprintf(myFile, "We begin with our original matrix:\n");
	texmatrix(myFile, A);

	/* Loop through columns in outer loop */
	for (int j = 0; j < n; j++) {
		
		/* Swap rows to position anchor */
		anchorExists = FALSE;
		for (int a = anchor; a < m; a++) {
			if (fabs(A[a][j]) > EPSILON) {

				swap(A, a, anchor);
				anchorExists = TRUE;

				if (anchor != a) {
					fprintf(myFile, "We will swap Row %d with Row %d as a suitable pivot:\n", a+1, anchor+1);
					texmatrix(myFile, A);
				}

				break;

			}
		}

		/* If anchor doesn't exist, move onto next column */
		if (anchorExists == FALSE) {
			fprintf(myFile, "We skip column %d because no pivot (i.e. nonzero entry) exists in this column.\n", j+1);
			continue;
		}

		/* Normalize anchor row */
		mult(A, anchor, 1.0 / A[anchor][j]);
		fprintf(myFile, "We now normalize Row %d so the pivot becomes equal to 1:\n", anchor+1);
		texmatrix(myFile, A);

		/* Loop through rows to scalar add by anchor i.e. "eliminate" */
		for (int i = 0; i < m; i++) {

			if (i != anchor) {

				fprintf(myFile, "We now add Row %d multiplied by a factor of %.2f to Row %d.", anchor+1, -A[i][j], i+1);
				fprintf(myFile, "This eliminates the entry in Row %d for Column %d.", i+1, j+1);

				add(A, i, anchor, -A[i][j]);

				texmatrix(myFile, A);

			}

		}

		/* Since anchor exists, inc anchor for next col */
		anchor++;

		/* If columns don't "run out" but rows do */
		if (anchor > m) {
			break;
		}

	}

	fprintf(myFile, "And thus we have our matrix in its RREF form:\n");
	texmatrix(myFile, A);

	endambles(myFile);
	
	fclose(myFile);

}

/*
 * Row operation: swaps two rows
 * Parameters: A (double[m][n]) - the matrix to be operated upon
 *             a (int)          - row # of one row to swap
 *             b (int)          - row # of the other row to swap
 */
void swap(double A[m][n], int a, int b) {

	/* Create a temporary array to store row a */
	double temp[n];

	for (int i = 0; i < n; i++) {

		temp[i]   = A[a][i]; // Copy row a onto temp
		A[a][i] = A[b][i]; // Copy row b onto a
		A[b][i] = temp[i]; // Copy temp onto b

	}

}

/*
 * Row operation: multiplies a row by a constant
 * Parameters: A (double[m][n]) - the matrix to be operated upon
 *             a (int)          - the row to be operated upon
 *             c (double)       - scalar multiple to multiply row a by
 */
void mult(double A[m][n], int a, double c) {

	/* Loop through columns of row a */
	for (int i = 0; i < n; i++) {

		A[a][i] *= c; // multiply element by c

	}

}

/*
 * Row operation: adds a constant multiple of a row onto another
 * Parameters: A (double[m][n]) - the matrix to be operated upon
 *             a (int)          - the row to be added onto (a = a + cb)
 *             b (int)          - the 'unmodified' row
 *             c (double)       - the scalar multiplier
 */
void add(double A[m][n], int a, int b, double c) {

	/* Loop through the columns */
	for (int i = 0; i < n; i++) {

		A[a][i] += c * A[b][i]; // add scalar multiple of b onto a

	}

}

/*
 * Writes the intro stuff to the TeX file.
 * Parameters: myFile (FILE*) - the file to be written to
 */
void preambles(FILE* myFile) {

	fprintf(myFile, "\\documentclass{article}\n");
	fprintf(myFile, "\\usepackage[utf8]{inputenc}\n");
	fprintf(myFile, "\\usepackage{amsmath}\n\n");

	fprintf(myFile, "\\title{Gaussian-Jordan Elimination of a $%d \\times %d$ Matrix}\n", m, n);
	fprintf(myFile, "\\author{%s}\n", TEXAUTHOR);
	fprintf(myFile, "\\date{%s}\n\n", TEXDATE);

	fprintf(myFile, "\\begin{document}\n\n");

	fprintf(myFile, "\\maketitle\n\n");

}

/*
* Writes the end stuff to the TeX file.
* Parameters: myFile (FILE*) - the file to be written to
*/
void endambles(FILE* myFile) {

	fprintf(myFile, "\\end{document}");

}

/*
 * Writes a given matrix in TeX to file
 * Parameters: myFile (FILE*)        - the file to be written to
 *             A      (double[m][n]) - the matrix to write
 */
void texmatrix(FILE* myFile, double A[m][n]) {

	fprintf(myFile, "\\[\n");
	fprintf(myFile, "\\begin{bmatrix}\n");

	for (int row = 0; row < m; row++) {

		/* Write initial entry */
		fprintf(myFile, "%.2f", A[row][0]);

		/* Loop through writing columns */
		for (int col = 1; col < n; col++) {
			fprintf(myFile, " & %.2f", A[row][col]);
		}

		/* Newline slashes */
		if (row != (m - 1)) {
			fprintf(myFile, " \\\\");
		}

		fprintf(myFile, "\n");

	}

	fprintf(myFile, "\\end{bmatrix}\n");
	fprintf(myFile, "\\]\n");

}