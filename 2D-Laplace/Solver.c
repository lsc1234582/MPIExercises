#include "Utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintHelp(void)
{
    printf("Usage: Solver <ParamsFile>\n");
}

int main(int argc, char**argv)
{
    printf("Info: Parsing args\n");
    /* Parse args */
    if (argc != 2)
    {
        PrintHelp();
        exit(1);
    }
    char paramsFileName[128];
    if(strlen(argv[1]) > 128)
    {
        printf("Error: Parameter file name too long\n");
        exit(1);
    }
    strncpy(paramsFileName, argv[1], strlen(argv[1]));
    Params params;
    if (ParseParameterFile(paramsFileName, &params))
    {
        printf("Error: Error in reading parameter file: %s\n", paramsFileName);
        exit(1);
    }
    printf("Info: Parameters:\n");
    PrintParameters(&params);

    /* Initialise grid by reading "initial.dat"
     * Set boundary values to those set by the analitycal solution (initial.dat), and all the interior values to zeros.
     */
    double** grid1 = malloc(params.m_NRow * sizeof(double*));
    for (int i = 0; i < params.m_NRow; ++i)
    {
        grid1[i] = calloc(params.m_NCol, sizeof(double));
    }
    double** grid2 = malloc(params.m_NRow * sizeof(double*));
    for (int i = 0; i < params.m_NRow; ++i)
    {
        grid2[i] = calloc(params.m_NCol, sizeof(double));
    }

    /* Parse initial.dat */
    printf("Info: Parsing initial.dat");
    FILE* fptr;
    if ((fptr = fopen("initial.dat", "r")) == NULL)
    {
        printf("Error: Cannot open initial.dat for reading\n");
        exit(1);
    }
    const double dx = (params.m_XMax - params.m_XMin) / (params.m_NRow - 1);
    const double dy = (params.m_YMax - params.m_YMin) / (params.m_NCol - 1);
    double x = params.m_XMin;
    double y = params.m_YMin;
    for (size_t i = 0; i < params.m_NRow; ++i)
    {
        y = params.m_YMin;
        for (size_t j = 0; j < params.m_NCol; ++j)
        {
            if(fscanf(fptr, "%*f %*f %lf\n", &grid1[i][j]) < 0)
            {
                exit(1);
            }
            grid2[i][j] = grid1[i][j];
            //printf("READ: %f %f %f\n", x, y, grid1[i][j]);
            y += dy;
        }
        if(fscanf(fptr, "\n") < 0)
        {
            exit(1);
        }
        x += dx;
    }
    if(fclose(fptr))
    {
        printf("Error: Cannot close initial.dat\n");
        exit(1);
    }

    /* Solve boundary value problem with Jacobi iteration method */
    printf("Info: Solving 2d Laplace with Jacobi iteration method\n");
    double maxDiff;
    do
    {
        maxDiff = 0.0;
        x = params.m_XMin;
        for (size_t i = 1; i < params.m_NRow - 1; ++i)
        {
            y = params.m_YMin;
            for (size_t j = 1; j < params.m_NCol - 1; ++j)
            {
                double term1 = (grid1[i-1][j] + grid1[i+1][j]) / (dx * dx) + (grid1[i][j-1] + grid1[i][j+1]) / (dy * dy); 
                double term2 = (dx * dx * dy * dy) / (2 * dx * dx + 2 * dy * dy);
                grid2[i][j] = term1 * term2;
                //printf("%f\t", grid2[i][j]);
                double diff = fabs(grid2[i][j] - grid1[i][j]);
                maxDiff = diff > maxDiff ? diff : maxDiff;
                y += dy;
            }
            //printf("\n");
            x += dx;
        }
        double** tempGrid = grid2;
        grid2 = grid1;
        grid1 = tempGrid;
        tempGrid = NULL;
        //printf("MAX_DIFF: %f\n", maxDiff);
    } while (maxDiff > params.m_Tolerance);

    /* Write results */
    printf("Info: Writing results to laplace.dat\n");
    FILE* wfptr;
    if ((wfptr = fopen("laplace.dat", "w")) == NULL)
    {
        printf("Error: Cannot open laplace.dat for writing\n");
        exit(1);
    }

    x = params.m_XMin;
    for (size_t i = 0; i < params.m_NRow; ++i)
    {
        y = params.m_YMin;
        for (size_t j = 0; j < params.m_NCol; ++j)
        {
            fprintf(wfptr, "%f %f %f\n", x, y, grid1[i][j]);
            y += dy;
        }
        fprintf(wfptr, "\n");
        x += dx;
    }
    if(fclose(wfptr))
    {
        printf("Error: Cannot close laplace.dat\n");
        exit(1);
    }

    /* Clean up */
    for (int i = 0; i < params.m_NRow; ++i)
    {
        free(grid2[i]);
    }
    free(grid2);
    for (int i = 0; i < params.m_NRow; ++i)
    {
        free(grid1[i]);
    }
    free(grid1);
    printf("Info: Exiting\n");
}
