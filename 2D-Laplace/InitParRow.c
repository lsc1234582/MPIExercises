#include "Utils.h"

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintHelp(void)
{
    printf("Usage: InitParRow <ParamsFile> <FunctionSelection>\n\
            <FunctionSelect>: [0..3]\n");
}

int main(int argc, char**argv)
{
    int size;
    int rank;
    int funcSelection;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Datatype ParamsMPIType;
    CreateGridParameterMPIStructDataType(&ParamsMPIType);
    GridParams params;
    /* Parse args */
    if (rank == MASTER_RANK)
    {
        printf("Info: Parsing args\n");
        if (argc != 3)
        {
            PrintHelp();
            exit(1);
        }
        char paramsFileName[MAX_FILE_NAME_LENGTH];
        if(strlen(argv[1]) > MAX_FILE_NAME_LENGTH)
        {
            printf("Error: Parameter file name too long\n");
            exit(1);
        }
        strncpy(paramsFileName, argv[1], strlen(argv[1]));
        char* endChr;
        funcSelection = strtol(argv[2], &endChr, 10);
        if (endChr == argv[2] || funcSelection > 3 || funcSelection < 0)
        {
            PrintHelp();
            exit(1);
        }
        /* Parse and broadcast parameters */
        if (ParseGridParameterFile(paramsFileName, &params))
        {
            printf("Error: Error in reading parameter file: %s\n", paramsFileName);
            exit(1);
        }
        printf("Info: Parameters:\n");
        PrintGridParameters(&params);
        printf("Info: Function selection: %d\n", funcSelection);
    }
    MPI_Bcast(&params, 1, ParamsMPIType, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(&funcSelection, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

    GridPatchParams horPatch;
    GetGridPatchParams(&params, size, rank, size, 1, &horPatch);

    double (*func)(double, double) = funcs[funcSelection];

    printf("Info: Writing initial bounadry values and analytical solutions to initial.MPI_%d.dat and solution.MPI_%d.dat\n",
            rank, rank);
    /* Produce both initial boundaries (writing to 'initial.dat'), and analytical solutions (writing to
     * 'solution.dat')) */
    double** gridInit = AllocateInitGridPatch(&horPatch);
    if (gridInit == NULL)
    {
        printf("Error: Error in allocating grid\n");
        exit(1);
    }
    double** gridSol = AllocateInitGridPatch(&horPatch);
    if (gridSol == NULL)
    {
        printf("Error: Error in allocating grid\n");
        exit(1);
    }
    double x = horPatch.m_PatchX;
    double y = horPatch.m_PatchY;
    for (size_t i = horPatch.m_AboveMargin; i < horPatch.m_NRow; ++i)
    {
        int globalI = horPatch.m_PatchI + i;
        y = horPatch.m_PatchY;
        for (size_t j = horPatch.m_LeftMargin; j < horPatch.m_NCol; ++j)
        {
            gridSol[i][j] = func(x, y);
            if (globalI > 0 && globalI < (params.m_NRow - 1) && j > 0 && j < (params.m_NCol - 1))
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
    sprintf(initialDatFileName, "initial.MPI_%d.dat", rank);
    char solutionDatFileName[MAX_FILE_NAME_LENGTH];
    sprintf(solutionDatFileName, "solution.MPI_%d.dat", rank);

    if (WriteGridPatch(initialDatFileName, &horPatch, gridInit))
    {
        exit(1);
    }
    if (WriteGridPatch(solutionDatFileName, &horPatch, gridSol))
    {
        exit(1);
    }
    /* Clean up */
    FreeGridPatch(&horPatch, gridInit);
    FreeGridPatch(&horPatch, gridSol);
    printf("Info: Exiting\n");
    MPI_Finalize();
}
