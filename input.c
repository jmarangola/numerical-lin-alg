#include <stdio.h>
#include <stdlib.h>

int read_square(size_t n, int **a, const char* file) {
    FILE *fp;
    fp = fopen(file, "r");
    if (fp == NULL) {
        printf("Input file not found.\n");
        return 0;
    }
    for(size_t i = 0; i < n; ++i) {
        for(size_t j = 0; j < n; ++j) {
            fscanf(fp, "%d", a[i] + j);
            printf("%d  ", a[i][j]);
        }
        printf("\n");
    }
    printf("\n");
    fclose (fp); 
    return 1; 
}

int main(void) {
    int **a = malloc(sizeof(int *) * 4);
    for (int i = 0; i < 4; i++) {
        a[i] = malloc(sizeof(int) * 4);
    }
    read_square(5, a, "input.txt");
}