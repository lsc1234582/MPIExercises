/* *
 * Solve laplace equation (parallel across rows and columns).
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
#define LEFT_TAG 3
#define RIGHT_TAG 4

void PrintHelp(void)
{
    printf("Usage: SolvePar <NumPatchInX> <NumPatchInY> <ParamsFile> [ResultDatBaseFileName [InitialDatBaseFileName]] \n");
}

int main(int argc, char**argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    int numPatchInX;
    int numPatchInY;
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
        if (argc < 4 || argc > 6)
        {
            PrintHelp();
            exit(1);
        }
        char* endChr;
        numPatchInX = strtol(argv[1], &endChr, 10);
        if (endChr == argv[1] || numPatchInX < 1)
        {
            PrintHelp();
            exit(1);
        }
        numPatchInY = strtol(argv[2], &endChr, 10);
        if (endChr == argv[2] || numPatchInY < 1)
        {
            PrintHelp();
            exit(1);
        }

        if(strlen(argv[3]) > MAX_FILE_NAME_LENGTH - 1)
        {
            printf("Error: Parameter file name too long\n");
            exit(1);
        }
        char paramsFileName[MAX_FILE_NAME_LENGTH];
        strncpy(paramsFileName, argv[3], MAX_FILE_NAME_LENGTH);

        if (argc > 4)
        {
            if(strlen(argv[4]) > MAX_FILE_NAME_LENGTH - 1)
            {
                printf("Error: Parameter file name too long\n");
                exit(1);
            }
            strncpy(resultDatBaseFileName, argv[4], MAX_FILE_NAME_LENGTH);
        }

        if (argc > 5)
        {
            if(strlen(argv[5]) > MAX_FILE_NAME_LENGTH - 1)
            {
                printf("Error: Parameter file name too long\n");
                exit(1);
            }
            strncpy(initialDatBaseFileName, argv[5], MAX_FILE_NAME_LENGTH);
        }

        if (ParseGridParameterFile(paramsFileName, &params))
        {
            exit(1);
        }

#ifdef DEBUG
        printf("Debug: Parameters:\n");
        PrintGridParameters(&params);
#endif
        printf("Info: Size: %d\n", size);
    }
    /* Broadcast common parameters to all processes */
    MPI_Bcast(&params, 1, ParamsMPIType, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(&numPatchInX, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(&numPatchInY, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(resultDatBaseFileName, MAX_FILE_NAME_LENGTH, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(initialDatBaseFileName, MAX_FILE_NAME_LENGTH, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);

    /* Allocate local grid patch */
    GridPatchParams patchParam;
    GetGridPatchParams(&params, size, rank, numPatchInX, numPatchInY, &patchParam);
    MPI_Datatype ColumnMarginElementT;
    CreateColumnMarginElementMPIDataType(&patchParam, &ColumnMarginElementT);
#ifdef DEBUG
    printf("Rank: %d, above_rank: %d, below_rank: %d\n", rank, patchParam.m_AboveRank, patchParam.m_BelowRank);
#endif
    if (patchParam.m_PatchI > 0)
    {
        assert(patchParam.m_AboveMargin == 1);
    }
    if (patchParam.m_PatchI < numPatchInX - 1)
    {
        assert(patchParam.m_BelowMargin == 1);
    }
    if (patchParam.m_PatchJ > 0)
    {
        assert(patchParam.m_LeftMargin == 1);
    }
    if (patchParam.m_PatchJ < numPatchInY - 1)
    {
        assert(patchParam.m_RightMargin == 1);
    }

    /* Read initial grid values */
    double** grid1 = AllocateInitGridPatch(&patchParam);
    if (grid1 == NULL)
    {
        pprintf("Error: Error in allocating grid\n");
        exit(1);
    }
    double** grid2 = AllocateInitGridPatch(&patchParam);
    if (grid2 == NULL)
    {
        pprintf("Error: Error in allocating grid\n");
        /* Clean up */
        FreeGrid(grid1);
        exit(1);
    }

    char initialDatFileName[MAX_FILE_NAME_LENGTH];
    snprintf(initialDatFileName, MAX_FILE_NAME_LENGTH, "%s.MPI_%d.dat", initialDatBaseFileName, rank);
    if (ReadGridPatch(initialDatFileName, &patchParam, grid1, grid2))
    {
        /* Clean up */
        FreeGrid(grid1);
        FreeGrid(grid2);
        exit(1);
    }

    /* Solve boundary value problem with Jacobi iteration method */
    pprintf("Info: Solving...\n");
    const double dx = patchParam.m_Dx;
    const double dy = patchParam.m_Dy;
    double maxDiff;
    double globalMaxDiff;
    double** tempGrid = NULL;
    int iterations = 0;
    do
    {
        /* Exchange boundary values */
        /* First exchange left and right columns horizontally */
        MPI_Request reqs[4];
        int numReqs = 0;
        MPI_Irecv(grid1[patchParam.m_AboveMargin], patchParam.m_NRow, ColumnMarginElementT, patchParam.m_LeftRank,
                RIGHT_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Isend(&grid1[patchParam.m_AboveMargin][patchParam.m_LeftMargin], patchParam.m_NRow, ColumnMarginElementT, patchParam.m_LeftRank,
                LEFT_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Irecv(&grid1[patchParam.m_AboveMargin][patchParam.m_NTotCol - 1], patchParam.m_NRow, ColumnMarginElementT, patchParam.m_RightRank,
                LEFT_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Isend(&grid1[patchParam.m_AboveMargin][patchParam.m_NTotCol - 2], patchParam.m_NRow, ColumnMarginElementT, patchParam.m_RightRank,
                RIGHT_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Waitall(numReqs, reqs, MPI_STATUSES_IGNORE);
        /* Then exchange top and bottom rows vertically, including any diagnal elements  */
        numReqs = 0;
        MPI_Irecv(grid1[0], patchParam.m_NTotCol, MPI_DOUBLE, patchParam.m_AboveRank,
                DOWN_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Isend(grid1[1], patchParam.m_NTotCol, MPI_DOUBLE, patchParam.m_AboveRank,
                UP_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Irecv(grid1[patchParam.m_NTotRow - 1], patchParam.m_NTotCol, MPI_DOUBLE, patchParam.m_BelowRank,
                UP_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Isend(grid1[patchParam.m_NTotRow - 2], patchParam.m_NTotCol, MPI_DOUBLE, patchParam.m_BelowRank,
                DOWN_TAG, MPI_COMM_WORLD, &reqs[numReqs++]);
        MPI_Waitall(numReqs, reqs, MPI_STATUSES_IGNORE);
        maxDiff = 0.0;
        for (size_t i = patchParam.m_AbovePadding; i < patchParam.m_NTotRow - patchParam.m_BelowPadding; ++i)
        {
            for (size_t j = patchParam.m_LeftPadding; j < patchParam.m_NTotCol - patchParam.m_RightPadding; ++j)
            {
                double term1 = (grid1[i-1][j] + grid1[i+1][j]) / (dx * dx) + (grid1[i][j-1] + grid1[i][j+1]) / (dy * dy);
                double term2 = (dx * dx * dy * dy) / (2 * dx * dx + 2 * dy * dy);
                grid2[i][j] = term1 * term2;
                double diff = fabs(grid2[i][j] - grid1[i][j]);
                maxDiff = diff > maxDiff ? diff : maxDiff;
            }
        }
        tempGrid = grid2;
        grid2 = grid1;
        grid1 = tempGrid;
        tempGrid = NULL;

        /* Get the global MaxDiff */
        MPI_Allreduce(&maxDiff, &globalMaxDiff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        //MPI_Barrier(MPI_COMM_WORLD);
        iterations++;
        //pprintf("MAX_DIFF: %f\n", maxDiff);
    } while (globalMaxDiff > params.m_Tolerance);

    /* Write results */
    char resultDatFileName[MAX_FILE_NAME_LENGTH];
    snprintf(resultDatFileName, MAX_FILE_NAME_LENGTH, "%s.MPI_%d.dat", resultDatBaseFileName, rank);
    WriteGridPatch(resultDatFileName, &patchParam, grid2);

    /* Clean up */
    FreeGrid(grid1);
    FreeGrid(grid2);

    pprintf("Info: Exiting\n");
    MPI_Finalize();
}
