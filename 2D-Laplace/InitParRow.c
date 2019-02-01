#include "Utils.h"

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void AllocateGridHorPatch(const Params* params, GridPatch* horPatch)
{
    int size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    const int fullRowSize = ceil(params->m_NRow / size);
    const int partialRowSize = params->m_NRow % size;
    horPatch->m_NRow = (partialRowSize == 0 || rank < size - 1) ? fullRowSize : partialRowSize;
    horPatch->m_PatchI = rank * fullRowSize;
    const double dx = GetDx(params);
    const double dy = GetDy(params);
    horPatch->m_PatchX = params->m_XMin + horPatch->m_PatchI * dx;
    horPatch->m_NCol = params->m_NCol;
    horPatch->m_PatchJ = 0;
    horPatch->m_PatchY = params->m_YMin;
    horPatch->m_Dx = dx;
    horPatch->m_Dy = dy;
}

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
    CreateParameterMPIStructDataType(&ParamsMPIType);
    Params params;
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
        if (ParseParameterFile(paramsFileName, &params))
        {
            printf("Error: Error in reading parameter file: %s\n", paramsFileName);
            exit(1);
        }
        printf("Info: Parameters:\n");
        PrintParameters(&params);
        printf("Info: Function selection: %d\n", funcSelection);
    }
    MPI_Bcast(&params, 1, ParamsMPIType, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(&funcSelection, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

    GridPatch horPatch;
    AllocateGridHorPatch(&params, &horPatch);

    const double dx = GetDx(&params);
    const double dy = GetDy(&params);

    double (*func)(double, double) = funcs[funcSelection];

    printf("Info: Writing initial bounadry values and analytical solutions to initial.MPI_%d.dat and solution.MPI_%d.dat\n",
            rank, rank);
    /* Produce both initial boundaries (writing to 'initial.dat'), and analytical solutions (writing to
     * 'solution.dat')) */
    double** gridInit = AllocateInitGrid(horPatch.m_NRow, params.m_NCol);
    if (gridInit == NULL)
    {
        printf("Error: Error in allocating grid\n");
        exit(1);
    }
    double** gridSol = AllocateInitGrid(horPatch.m_NRow, params.m_NCol);
    if (gridSol == NULL)
    {
        printf("Error: Error in allocating grid\n");
        exit(1);
    }
    double x = horPatch.m_PatchX;
    double y = params.m_YMin;
    for (size_t i = 0; i < horPatch.m_NRow; ++i)
    {
        int globalI = horPatch.m_PatchI + i;
        y = params.m_YMin;
        for (size_t j = 0; j < params.m_NCol; ++j)
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
            y += dy;
        }
        x += dx;
    }

    char initialDatFileName[MAX_FILE_NAME_LENGTH];
    sprintf(initialDatFileName, "initial.MPI_%d.dat", rank);
    char solutionDatFileName[MAX_FILE_NAME_LENGTH];
    sprintf(solutionDatFileName, "solution.MPI_%d.dat", rank);

    if (WriteGridHorPatch(initialDatFileName, &params, &horPatch, gridInit))
    {
        exit(1);
    }
    if (WriteGridHorPatch(solutionDatFileName, &params, &horPatch, gridSol))
    {
        exit(1);
    }
    /* Clean up */
    FreeGrid(horPatch.m_NRow, gridInit);
    FreeGrid(horPatch.m_NRow, gridSol);
    printf("Info: Exiting\n");
    MPI_Finalize();
}
