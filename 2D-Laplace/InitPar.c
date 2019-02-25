#include "Utils.h"

#include <mpi.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintHelp(void)
{
    printf("Usage: InitPar <NumPatchInX> <NumPatchInY> <ParamsFile> <FunctionSelection> [InitDatBaseFileName [SolDatBaseFileName]]\n\
            <FunctionSelect>: [0..3]\n");
}

int main(int argc, char**argv)
{
    int size;
    int rank;
    int funcSelection;
    int numPatchInX, numPatchInY;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Datatype ParamsMPIType;
    CreateGridParameterMPIStructDataType(&ParamsMPIType);
    GridParams params;
    char initialDatBaseFileName[MAX_FILE_NAME_LENGTH] = "initial";
    char solutionDatBaseFileName[MAX_FILE_NAME_LENGTH] = "solution";
    /* Parse args */
    if (rank == MASTER_RANK)
    {
        printf("Info: Parsing args\n");
        if (argc < 5 || argc > 7)
        {
            PrintHelp();
            exit(1);
        }
        char paramsFileName[MAX_FILE_NAME_LENGTH];
        if(strlen(argv[3]) > MAX_FILE_NAME_LENGTH - 1)
        {
            printf("Error: Parameter file name too long\n");
            exit(1);
        }
        strncpy(paramsFileName, argv[3], MAX_FILE_NAME_LENGTH);
        char* endChr;
        funcSelection = strtol(argv[4], &endChr, 10);
        if (endChr == argv[4] || funcSelection > 3 || funcSelection < 0)
        {
            PrintHelp();
            exit(1);
        }
        numPatchInX = strtol(argv[1], &endChr, 10);
        if (endChr == argv[1])
        {
            PrintHelp();
            exit(1);
        }
        numPatchInY = strtol(argv[2], &endChr, 10);
        if (endChr == argv[2])
        {
            PrintHelp();
            exit(1);
        }
        assert(size == numPatchInX * numPatchInY);
        /* Parse and broadcast parameters */
        if (ParseGridParameterFile(paramsFileName, &params))
        {
            printf("Error: Error in reading parameter file: %s\n", paramsFileName);
            exit(1);
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

        if (argc > 6)
        {
            if(strlen(argv[6]) > MAX_FILE_NAME_LENGTH - 1)
            {
                printf("Error: Parameter file name too long\n");
                exit(1);
            }
            strncpy(solutionDatBaseFileName, argv[6], MAX_FILE_NAME_LENGTH);
        }
        printf("Info: Parameters:\n");
        PrintGridParameters(&params);
        printf("Info: Function selection: %d\n", funcSelection);
    }
    MPI_Bcast(&params, 1, ParamsMPIType, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(&funcSelection, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(&numPatchInX, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(&numPatchInY, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(initialDatBaseFileName, MAX_FILE_NAME_LENGTH, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(solutionDatBaseFileName, MAX_FILE_NAME_LENGTH, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);

    GridPatchParams horPatch;
    GetGridPatchParams(&params, size, rank, numPatchInX, numPatchInY, &horPatch);

    double (*func)(double, double) = funcs[funcSelection];

    /* Produce both initial boundaries, and analytical solutions */
    double** gridInit = AllocateInitGridPatch(&horPatch);
    if (gridInit == NULL)
    {
        printf("Error: Error in allocating grid\n");
        exit(1);
    }
    double** gridSol = AllocateInitGridPatch(&horPatch);
    if (gridSol == NULL)
    {
        /* Clean up */
        FreeGrid(gridInit);
        printf("Error: Error in allocating grid\n");
        exit(1);
    }
    double x = horPatch.m_PatchX;
    double y = horPatch.m_PatchY;
    for (size_t i = horPatch.m_AboveMargin; i < horPatch.m_NTotRow - horPatch.m_BelowMargin; ++i)
    {
        int globalI = horPatch.m_PatchI + i - horPatch.m_AboveMargin;
        y = horPatch.m_PatchY;
        for (size_t j = horPatch.m_LeftMargin; j < horPatch.m_NTotCol - horPatch.m_RightMargin; ++j)
        {
            int globalJ = horPatch.m_PatchJ + j - horPatch.m_LeftMargin;
            gridSol[i][j] = func(x, y);
            if (globalI > 0 && globalI < (params.m_NRow - 1) && globalJ > 0 && globalJ < (params.m_NCol - 1))
            {
                gridInit[i][j] = 0.0;
            }
            else
            {
                gridInit[i][j] = gridSol[i][j];
            }
            y += horPatch.m_Dy;
        }
        x += horPatch.m_Dx;
    }

    char initialDatFileName[MAX_FILE_NAME_LENGTH];
    snprintf(initialDatFileName, MAX_FILE_NAME_LENGTH, "%s.MPI_%d.dat", initialDatBaseFileName, rank);
    char solutionDatFileName[MAX_FILE_NAME_LENGTH];
    snprintf(solutionDatFileName, MAX_FILE_NAME_LENGTH, "%s.MPI_%d.dat", solutionDatBaseFileName, rank);

    WriteGridPatch(initialDatFileName, &horPatch, gridInit);
    WriteGridPatch(solutionDatFileName, &horPatch, gridSol);

    /* Clean up */
    FreeGrid(gridInit);
    FreeGrid(gridSol);
    printf("Info: Exiting\n");
    MPI_Finalize();
}
