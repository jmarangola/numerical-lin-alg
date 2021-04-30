#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int read_square(size_t n, int **a, const char *file) {
    FILE *fp;
    fp = fopen(file, "r");
    if (fp == NULL) {
        printf("Input file not found.\n");
        return 0;
    }
    for(size_t i = 0; i < n; ++i) {
        for(size_t j = 0; j < n; ++j) {
            fscanf(fp, "%d", a[i] + j);
        }
        printf("\n");
    }
    printf("\n");
    fclose (fp); 
    return 1; 
}

int write_square(size_t n, int **a, const char *file) {
    FILE *fp;
    fp = fopen(file, "w");
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            fprintf(fp, "%d ", a[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 1;
}


int main(void) {
    int **a = malloc(sizeof(int *) * 4);
    for (int i = 0; i < 4; i++) {
        a[i] = malloc(sizeof(int) * 4);
    }
    double *d = malloc(sizeof(double) * 10);
    read_arbitrary("input.txt", d);
    write_square(4, a, "output.txt");
    
}