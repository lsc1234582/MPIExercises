#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"

#define MAX_VALUE 10

void print_help(void)
{
    printf("Format: main <Number of Row> <Number of Col>\n");
}

int main(int argc, char** argv)
{
    if(argc < 3)
    {
        print_help();
        exit(1);
    }

    const int AROW = atoi(argv[1]);
    const int ACOL = atoi(argv[2]);

    int** a = (int**) malloc(AROW * ACOL * sizeof(int *));
    for (int i = 0; i < AROW; ++i)
    {
        a[i] = (int*) calloc(ACOL, sizeof(int));
    }
    int* b = (int*) calloc(ACOL, sizeof(int));
    int* c = (int*) calloc(AROW, sizeof(int));

    /* Generating Random Values for A & B Array*/
    srand(time(NULL));
    for (int i=0;i<AROW;i++)
    {
        for (int j=0;j<ACOL;j++)
        {
            if (i==0) b[j] = rand() % MAX_VALUE;
            a[i][j] = rand() % MAX_VALUE;
        }
    }

    /* Printing the Matrix*/
#ifdef DEBUG
    printf("Matrix A :\n");
    for (int i=0;i<AROW;i++)
    {
        for (int j=0;j<ACOL;j++)
        {
            printf("%3d ", a[i][j]);
        }
        printf("\n");
    }
    printf("nMatrix B :\n");
    for (int i=0;i<ACOL;i++)
    {
        printf("%3d ", b[i]);
    }
    printf("\n\n");
#endif

    /* (3) Gathering the result from other processes*/
    for (int i=0;i<AROW;i++)
    {
        for (int j=0; j < ACOL; ++j)
        {
            c[i] += a[i][j] * b[j];
        }
#ifdef DEBUG
        printf("P%d : c[%d]\t= %d\n", rank, i, c[i]);
#endif
    }

    for (int i = 0; i < AROW; ++i)
    {
        free(a[i]);
    }
    free(a);
    free(b);
    free(c);
    return 0;
}
