/* *
 * Solve laplace equation (parallel across rows).
 */
#include "Utils.h"

#include <mpi.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define UP_TAG 1
#define DOWN_TAG 2

void PrintHelp(void)
{
    printf("Usage: SolveParRow <ParamsFile> [ResultDatBaseFileName [InitialDatBaseFileName]] \n");
}

int main(int argc, char**argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Datatype ParamsMPIType;
    CreateGridParameterMPIStructDataType(&ParamsMPIType);
    GridParams params;
    char resultDatBaseFileName[MAX_FILE_NAME_LENGTH] = "laplace";
    char initialDatBaseFileName[MAX_FILE_NAME_LENGTH] = "initial";

    pprintf("Info: Command ");
    PrintCLArgs(argc, argv);
    if (rank == MASTER_RANK)
    {
        /* Parse args */
        printf("Info: Parsing args\n");
        if (argc < 2 || argc > 4)
        {
            PrintHelp();
            exit(1);
        }
        char paramsFileName[MAX_FILE_NAME_LENGTH];
        if(strlen(argv[1]) > MAX_FILE_NAME_LENGTH - 1)
        {
            printf("Error: Parameter file name too long\n");
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


        if (ParseGridParameterFile(paramsFileName, &params))
        {
            exit(1);
        }
        printf("Info: Parameters:\n");
        PrintGridParameters(&params);
        printf("Info: Size: %d\n", size);
    }
    /* Broadcast common parameters to all processes */
    MPI_Bcast(&params, 1, ParamsMPIType, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(resultDatBaseFileName, MAX_FILE_NAME_LENGTH, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(initialDatBaseFileName, MAX_FILE_NAME_LENGTH, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);

    /* Allocate local grid patch */
    GridPatchParams horPatch;
    GetGridPatchParams(&params, size, rank, size, 1, &horPatch);
    assert(horPatch.m_LeftMargin == 0 && horPatch.m_RightMargin == 0);
    printf("Rank: %d, above_rank: %d, below_rank: %d\n", rank, horPatch.m_AboveRank, horPatch.m_BelowRank);
    if (horPatch.m_PatchI > 0)
    {
        printf("AbovemARGIN: %d %d %d\n", horPatch.m_AboveMargin, size, horPatch.m_PatchI);
        assert(horPatch.m_AboveMargin == 1);
    }
    if (horPatch.m_PatchI < size - 1)
    {
        printf("BelowmARGIN: %d", horPatch.m_BelowMargin);
        assert(horPatch.m_BelowMargin == 1);
    }

    /* Read initial grid values */
    double** grid1 = AllocateInitGridPatch(&horPatch);
    if (grid1 == NULL)
    {
        pprintf("Error: Error in allocating grid\n");
        exit(1);
    }
    double** grid2 = AllocateInitGridPatch(&horPatch);
    if (grid2 == NULL)
    {
        pprintf("Error: Error in allocating grid\n");
        /* Clean up */
        FreeGrid(grid1);
        exit(1);
    }

    char initialDatFileName[MAX_FILE_NAME_LENGTH];
    snprintf(initialDatFileName, MAX_FILE_NAME_LENGTH, "%s.MPI_%d.dat", initialDatBaseFileName, rank);
    if (ReadGridPatch(initialDatFileName, &horPatch, grid1, grid2))
    {
        /* Clean up */
        FreeGrid(grid1);
        FreeGrid(grid2);
        exit(1);
    }

    /* Solve boundary value problem with Jacobi iteration method */
    pprintf("Info: Solving...\n");

    const double dx = horPatch.m_Dx;
    const double dy = horPatch.m_Dy;
    double x = horPatch.m_PatchX;
    double y = horPatch.m_PatchY;
    double maxDiff;
    double globalMaxDiff;
    double** tempGrid = NULL;
    int iterations = 0;
    do
    {
        /* Exchange boundary values */
        MPI_Request reqs[4];
        int numReqs = 0;
        MPI_Irecv(grid1[0], horPatch.m_NTotCol, MPI_DOUBLE, horPatch.m_AboveRank,
                DOWN_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Isend(grid1[1], horPatch.m_NTotCol, MPI_DOUBLE, horPatch.m_AboveRank,
                UP_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Irecv(grid1[horPatch.m_NTotRow - 1], horPatch.m_NTotCol, MPI_DOUBLE, horPatch.m_BelowRank,
                UP_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Isend(grid1[horPatch.m_NTotRow - 2], horPatch.m_NTotCol, MPI_DOUBLE, horPatch.m_BelowRank,
                DOWN_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Waitall(numReqs, reqs, MPI_STATUSES_IGNORE);
        maxDiff = 0.0;
        x = horPatch.m_PatchX;
        for (size_t i = horPatch.m_AbovePadding; i < horPatch.m_NTotRow - horPatch.m_BelowPadding; ++i)
        {
            y = horPatch.m_PatchY;
            for (size_t j = horPatch.m_LeftPadding; j < horPatch.m_NTotCol - horPatch.m_RightPadding; ++j)
            {
                double term1 = (grid1[i-1][j] + grid1[i+1][j]) / (dx * dx) + (grid1[i][j-1] + grid1[i][j+1]) / (dy * dy);
                double term2 = (dx * dx * dy * dy) / (2 * dx * dx + 2 * dy * dy);
                grid2[i][j] = term1 * term2;
                //pprintf("%f\t", grid2[i][j]);
                double diff = fabs(grid2[i][j] - grid1[i][j]);
                maxDiff = diff > maxDiff ? diff : maxDiff;
                y += dy;
            }
            //pprintf("\n");
            x += dx;
        }
        tempGrid = grid2;
        grid2 = grid1;
        grid1 = tempGrid;
        tempGrid = NULL;

        /* Get the global MaxDiff */
        MPI_Allreduce(&maxDiff, &globalMaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        /*
        if (iterations % 100 == 0)
        {
            char resultDatFileName[MAX_FILE_NAME_LENGTH];
            sprintf(resultDatFileName, "laplace.MPI_%d_iter_%d.dat", rank, iterations);
            if (WriteGridPatch(resultDatFileName, &horPatch, grid2))
            {
                exit(1);
            }
        }
        */
        //MPI_Barrier(MPI_COMM_WORLD);
        iterations++;
        //pprintf("MAX_DIFF: %f\n", maxDiff);
    } while (globalMaxDiff > params.m_Tolerance);

    /* Write results */
    char resultDatFileName[MAX_FILE_NAME_LENGTH];
    snprintf(resultDatFileName, MAX_FILE_NAME_LENGTH, "%s.MPI_%d.dat", resultDatBaseFileName, rank);
    WriteGridPatch(resultDatFileName, &horPatch, grid2);

    /* Clean up */
    FreeGrid(grid1);
    FreeGrid(grid2);

    pprintf("Info: Exiting\n");
    MPI_Finalize();
}
