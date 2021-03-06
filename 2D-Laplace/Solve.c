/* *
 * Solve laplace equation serially.
 */
#include "Utils.h"

#include <assert.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintHelp(void)
{
    pprintf("Usage: Solve <ParamsFile> [ResultDatBaseFileName [InitialDatBaseFileName]] \n");
}

int main(int argc, char**argv)
{
    char resultDatBaseFileName[MAX_FILE_NAME_LENGTH] = "laplace";
    char initialDatBaseFileName[MAX_FILE_NAME_LENGTH] = "initial";

    /* Parse args */
    pprintf("Info: Command ");
    PrintCLArgs(argc, argv);
    pprintf("Info: Parsing args\n");
    if (argc < 2 || argc > 4)
    {
        PrintHelp();
        exit(1);
    }
    char paramsFileName[MAX_FILE_NAME_LENGTH];
    if(strlen(argv[1]) > MAX_FILE_NAME_LENGTH - 1)
    {
        pprintf("Error: Parameter file name too long\n");
        exit(1);
    }
    strncpy(paramsFileName, argv[1], MAX_FILE_NAME_LENGTH);
    if (argc > 2)
    {
        if(strlen(argv[2]) > MAX_FILE_NAME_LENGTH - 1)
        {
            printf("Error: Parameter file name too long\n");
            exit(1);
        }
        strncpy(resultDatBaseFileName, argv[2], MAX_FILE_NAME_LENGTH);
    }

    if (argc > 3)
    {
        if(strlen(argv[3]) > MAX_FILE_NAME_LENGTH - 1)
        {
            printf("Error: Parameter file name too long\n");
            exit(1);
        }
        strncpy(initialDatBaseFileName, argv[3], MAX_FILE_NAME_LENGTH);
    }

    GridParams params;
    InitialiseGridParams(&params);
    if (ParseGridParameterFile(paramsFileName, &params))
    {
        pprintf("Error: Error in reading parameter file: %s\n", paramsFileName);
        exit(1);
    }
    assert(params.m_NCol > 1);
    assert(params.m_NRow > 1);
#ifdef DEBUG
    pprintf("Debug: Parameters:\n");
    PrintGridParameters(&params);
#endif

    /* Read initial grid values */
    double** grid1 = AllocateInitGrid(params.m_NRow, params.m_NCol);
    if (grid1 == NULL)
    {
        pprintf("Error: Error in allocating grid\n");
        exit(1);
    }
    double** grid2 = AllocateInitGrid(params.m_NRow, params.m_NCol);
    if (grid2 == NULL)
    {
        pprintf("Error: Error in allocating grid\n");
        /* Clean up */
        FreeGrid(grid1);
        exit(1);
    }

    char initialDatFileName[MAX_FILE_NAME_LENGTH];
    snprintf(initialDatFileName, MAX_FILE_NAME_LENGTH, "%s.dat", initialDatBaseFileName);
    if (ReadGrid(initialDatFileName, &params, grid1, grid2))
    {
        /* Clean up */
        FreeGrid(grid1);
        FreeGrid(grid2);
        exit(1);
    }

    /* Solve boundary value problem with Jacobi iteration method */
    pprintf("Info: Solving...\n");
    const double dx = GetDx(&params);
    const double dy = GetDy(&params);
    double** diffs = AllocateInitGrid(params.m_NRow, params.m_NCol);
    double maxDiff = params.m_Tolerance;
    int diffsCount = 0;
    double** tempGrid = NULL;
    do
    {
        diffsCount = 0;
        for (size_t i = 1; i < params.m_NRow - 1; ++i)
        {
            for (size_t j = 1; j < params.m_NCol - 1; ++j)
            {
                double term1 = (grid1[i-1][j] + grid1[i+1][j]) / (dx * dx) + (grid1[i][j-1] + grid1[i][j+1]) / (dy * dy);
                double term2 = (dx * dx * dy * dy) / (2 * dx * dx + 2 * dy * dy);
                grid2[i][j] = term1 * term2;
                diffs[0][diffsCount++] = fabs(grid2[i][j] - grid1[i][j]);
            }
        }

        maxDiff = params.m_Tolerance;
        for (size_t i = 0; i < diffsCount; ++i)
        {
            maxDiff = diffs[0][i] > maxDiff ? diffs[0][i] : maxDiff;
        }

        tempGrid = grid2;
        grid2 = grid1;
        grid1 = tempGrid;
        tempGrid = NULL;
    } while (maxDiff > params.m_Tolerance);

    /* Write results */
    char resultDatFileName[MAX_FILE_NAME_LENGTH];
    snprintf(resultDatFileName, MAX_FILE_NAME_LENGTH, "%s.dat", resultDatBaseFileName);
    WriteGrid(resultDatFileName, &params, grid2);

    /* Clean up */
    FreeGrid(grid1);
    FreeGrid(grid2);

    pprintf("Info: Exiting\n");
}
