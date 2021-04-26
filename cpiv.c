/*
    * John Marangola    
    * Inverse, Determinant, Rank
*/

#include <math.h>     
#include <stdio.h>                     
#include <stdlib.h>


void print_matrix(double **, int);
void print_col_vector(double *, int);

enum output_type {basic, verbose};

// Parameters for program
const int M = 5;
const int N = 5;
const double eps = 1e-10;
const enum output_type output = verbose;

void free_matrix(double **mat, int cols) {
    for (int i = 0; i < cols; i++)
        free(mat[i]);
    free(mat);
}
/**
 * @brief Compute the LU decomposition of an nxn square, real matrix a in place using pivoting.
 * The pivoting strategy is preserved within the unitary permutation matrix *p, which has an (n+1)th
 * element equal to the number of row/col exchanges performed in the decomposition.
 * @param epsilon_zero Minimum threshold for non-singularity based on numerical accumulation 
 * @c O(n^3)
 * @return <int> returns 1 if a is succesfully decompose, otherwise 0.
 */
int lu_decomposition_with_pivoting(double **a, int n, double epsilon_zero, int *p) {
    int k, i, j, r_max_index;
    double r_max, *ptr, absval;
    // Initialize permutation matrix
    for (i = 0; i <= n; i++)
        p[i] = i;
    // Begin LU Decomposition with pivoting
    for (i = 0; i < n; i++) {
        r_max = 0.0;
        r_max_index = i;
        for (k = i; k < n; k++) {
            if ((absval = fabs(a[k][i])) > r_max) {
                r_max = absval;
                r_max_index = k;
            }
        }
        // Check if singular within numerical bounds of epsilon_zero
        if (r_max < epsilon_zero) {
            printf("Matrix is singular.\n");
            return 0; 
        }
        if (r_max_index != i) {
            // Update permutation "psuedo"-matrix P
            j = p[i];
            p[i] = p[r_max_index];
            p[r_max_index] = j;
            // Perform pivoting in place on a
            ptr = a[i];
            a[i] = a[r_max_index];
            a[r_max_index] = ptr;
            // Update pivot count (n+1)th element to use later for determinant sign
            p[n]++;
        }
        // For j > i above the diagonal 
        for (j = i + 1; j < n; j++) { 
            a[j][i] /= a[i][i];
            for (k = i + 1; k < n; k++)
                a[j][k] -= a[j][i] * a[i][k];
        }
    }
    return 1;
}

/**
 * @brief LU decomposition without pivoting, this method is numerically unstable and should be avoided if possible.
 * Because of its instability, it will not always converge and numerical accumulation errors effect accuracy.
 **/ 
void lu_decomposition_without_pivoting(double **a, int n) {
    int i, j, k;
    // Crout's method without partial pivoting
    for(j = 0; j < n; j++) {
        for(i = 0 ; i <= j; i++) {
            if(i > 0) {
                for(k = 0; k < i; k++)
                    a[i][j] -= a[i][k] * a[k][j];
            }
        }
        if(j < n - 1) {
            for(i = j + 1; i < n; i++) {
                if(j > 0) {
                    for(k = 0; k < j; k++)
                        a[i][j] -= a[i][k] * a[k][j];
                }     
                a[i][j] = a[i][j] / a[j][j];
            }
        }
    }
}

/**
 * @brief Separate lower and upper triangular matrices l and u, respectively 
 * from a mutated matrix where lu decomp was performed "in place". Merely 
 * implemented for the purposes of visualizing solutions.
 */
void lower_upper(double **a, double **lower, double **upper, int n) {
    int i, j;
    for (i =0 ; i < n; i++) {
        for (j = 0; j < i; j++) {
            lower[i][j] = a[i][j];
        }
        lower[j][j] = 1; 
        for (j = i; j < n; j++) 
            upper[i][j] = a[i][j];
    }
}

/**
 * @brief Computes and stores the inverse of a real, square (nxn) matrix a in
 * inverse. Utilizes LU decomposition with pivoting, forward substitution and back substitution
 * to solve for the inverse of a matrix providing that it is nonsingular.
 * @c O(n^3)
 * @return <int> 1 if a is invertible, 0 otherwise.
 */
int invert(double **a, double **inverse, int n) {
    int *perm = (int *)malloc(sizeof(int) * n);
    lu_decomposition_with_pivoting(a, n, 1e-20, perm);
    // j = perm[i] is row of original matrix in the ith row of the matrix
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            if (perm[i] == j)
                inverse[i][j] = 1.0;
            else
                inverse[i][j] = 0.0;
            // Forward Substitution Lx = b_i
            for (int k = 0; k < i; k++)
                inverse[i][j] -= a[i][k] * inverse[k][j];
        }
        // Back substitution Ux = b_i
        for (int i = n - 1; i >= 0; i--) {
            for (int k = i + 1; k < n; k++)
                inverse[i][j] -= a[i][k] * inverse[k][j];
            inverse[i][j] = inverse[i][j] / a[i][i];
        }
    }
    return 1;
}

// Compute the sign of a scalar
int sign(double d) { return ((d >= 0) ? 1: -1); }

/**
 * @brief Sets id to the nxn identity
 */
void eye(double **id, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) 
                id[i][j] = 1;
            else 
                id[i][j] = 0;
        }
    }
}

/**
 * @brief Compute squared norm of a real-valued n-dimensional vector.
 * @return <double>
 */
double norm_squared(double *x, int n) {
    double value = 0.0;
    for (int i = 0; i < n; i++) value += x[i]*x[i];
    return value;
}

/**
 * @brief Set all of the values of an n-dimensional vector x to zero
 */
void reset(double *x, int n) {
    for (int i = 0; i < n; i++) x[i] = 0.0;
}

/**
 * @brief Prints an mxn matrix a to four decimal places of precision
 */
void print_rect(double **a, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%.4f  ", a[i][j]);
        }
        putchar('\n');
    }
    putchar('\n');
}

/**
 * @brief Computes the rank of a real-valued mxn matrix a.
 * Method utilizes a rank-revealing QR factorization (RRQR)
 * @param epsilon tolerance to numerical accumulation 
 * @return <int> rank of matrix  
 */
int rank(double **a, int m, int n, double epsilon) {
    double z[m], v[m], sum_sqs, sum;
    for (int k = 0; k < m; k++) {
        reset(v, m);
        reset(z, m);
        sum_sqs = 0.0;
        // Slice z from current col of matrix a (z = A_{:k})
        for (int i = k; i < m; i++) 
            z[i - k] = a[i][k];
        // Determine the Householder reflector
        v[0] = sign(z[0]) * sqrt(norm_squared(z, m - k)) + z[0];
        sum_sqs += v[0] * v[0];
        for (int i = 1; i < m - k; i++) {
            v[i] = z[i];
            sum_sqs += v[i]*v[i];
        }
        for (int i = 0; i < m; i++) 
            v[i] = v[i]/sqrt(sum_sqs);
        // Apply the Householder transformation to the current matrix
        for (int j = k; j < m; j++) {
            sum = 0.0;
            for (int i = k; i < n; i++)
                sum += v[i - k] * a[i][j];
            for (int i = k; i < n; i++) {
                a[i][j] -= 2 * v[i - k] * sum;
            }
            if (output == verbose) print_rect(a, m, n);
        }    
    }
    // Determine rank of matrix with given tolerance epsilon
    int rank = 0;
    for (int i = 0; i < n; i++)
        if (fabs(a[i][i]) > epsilon) rank++;
    return rank;
}

/**
 * @brief Computes the determinant for a matrix A, U, L ∈ R^{nxn}
 * Provided with the LU decomposition of a matrix, the computation 
 * of its determinant is O(n) because det(A) = det(LU) = det(L)det(U) = det(U)
 * since all the diagonal elements of L are taken to be unity and 
 * the determinant of a triangular matrix is O(n).
 * @return <double> determinant if nonsingular, zero otherwise
 */
double det(double **a, int n) {
    int *perm = (int *)malloc(sizeof(int) * n);
    lu_decomposition_with_pivoting(a, n, 1e-20, perm);
    double det = a[0][0];
    int tperm = perm[n] - n;
    for (int i = 1; i < n; i++)
        det *= a[i][i];
    free(perm);
    // Compute determinant depending on swap parity:
    return (!(tperm % 2)) ? det : -det;
}

/**
 * @brief Outputs an n-dim column vector x
 */
void print_col_vector(double *x, int n) {
    for (int i = 0; i < n; i++) 
        printf("%f\n", x[i]);
    putchar('\n');
}

/**
 * @brief Prints an nxn real matrix to four decimal places
 */
void print_matrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j  < n; j++) {
            printf("%0.4f  ", matrix[i][j]);
        }
        putchar('\n');
    }
    putchar('\n');
}

/**
 * @brief Perform standard matrix multiplication a times b storing result in c.
 * a, b, c are all assumed to be nxn, real, (square) matrices.
 */
void multiply(double **a, double **b, double **c, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            c[i][j] = 0.0;
            for (int k = 0; k < n; k++)
                c[i][j] += a[i][k]*b[k][j];
        }
    }
}

int main(void) {

    

    // Initialize MxN real matrix rmat
    double rect[M][N] = {
        {8, 1, 6, 2, 1}, 
        {3, 5, 7, 5, 6}, 
        {4, 9, 2, 1, 1},
        {2, 1, 6, 1, 6}, 
        {2, 2, 2, 2, 1}
    };

    // Initialize NxN square matrix 
    double data[][N] = {
        {1.1, 4.0, 9, 2, 16},
        {17.2, -2.1, 4.2, 1.8, -6},
        {4.9, 6, -5.1, 2, 17.8},
        {-2.5, 11, 13, -14.2, 0.2},
        {1, 0, 2, 1, 0}
    };

    // Allocate memory
    double **rmat = (double **)malloc(sizeof(double *) * N);
    double **matrix = (double **) malloc(sizeof(double *) * N);
    double **lower = (double **) malloc(sizeof(double *) * N);
    double **upper = (double **) malloc(sizeof(double *) * N);
    double **product = (double **) malloc(sizeof(double *) * N);
    double **product2 = (double **) malloc(sizeof(double *) * N);
    double **copy1 = (double **)malloc(sizeof(double *) * N);
    double **copy2 = (double **)malloc(sizeof(double *) * N);
    double **inverse = (double **)malloc(sizeof(double *) * N);
    double **abs_error = (double **)malloc(sizeof(double *) * N);
    for (int i = 0; i < M; i++) 
        rmat[i] = rect[i];
    for (int i = 0; i < N; i++) {
        // Copy initialized matrix above
        matrix[i] = data[i];
        // Allocate memory for empty matrices
        lower[i] = (double *)malloc(sizeof(double) * N);
        upper[i] = (double *)malloc(sizeof(double) * N);
        product[i] = (double *)malloc(sizeof(double) * N);
        product2[i] = (double *)malloc(sizeof(double) * N);
        inverse[i] = (double *)malloc(sizeof(double) * N);
        copy1[i] = (double *)malloc(sizeof(double) * N);
        copy2[i] = (double *)malloc(sizeof(double) * N);
        abs_error[i] = (double *)malloc(sizeof(double) * N);

    }

    // Copy matrix to copy
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++){
            copy1[i][j] = matrix[i][j];
            copy2[i][j] = matrix[i][j];
        }
    }

    int *P = malloc(sizeof(int) * N); // Unitary permutation matrix 
    printf("matrix: \n");
    print_matrix(matrix, N);
    printf("rmat \n");
    print_rect(rmat, M, N);
    
    
    if (output == verbose) {
        printf("Decomposing matrix into LU factorization: \n");
        lu_decomposition_with_pivoting(copy1, N, 1e-20, P);
        lower_upper(copy1, lower, upper, N);
        printf("L:\n");
        print_matrix(lower, N);
        printf("U:\n");
        print_matrix(upper, N);
        multiply(lower, upper, product, N);
        printf("L times U: \n");
        print_matrix(product, N);
    }
    invert(copy2, inverse, N);
    printf("Inverse: \n");
    print_matrix(inverse, N);
    if (output == verbose) {
        printf("A^-1*A = I:\n");
        multiply(inverse, matrix, product2, N);
        print_matrix(product2, N);
    }
    printf("det(matrix) = %f\n", det(matrix, N));
    printf("Rank(rmat): %i\n", rank(matrix, M, N, eps));
   
    free_matrix(lower, N);
    free(upper);
    free(matrix);
    free(product);
    free(copy2);
    free(P);
    free(inverse);
}