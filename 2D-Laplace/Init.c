/* *
 * Initialise boundary values serially.
 */
#include "Utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintHelp(void)
{
    printf("Usage: Init <ParamsFile> <FunctionSelection> [InitDatBaseFileName [SolDatBaseFileName]]\n\
            <FunctionSelect>: [0..3]\n");
}

int main(int argc, char**argv)
{
    /* Parse args */
    pprintf("Info: Command ");
    PrintCLArgs(argc, argv);
    char initialDatBaseFileName[MAX_FILE_NAME_LENGTH] = "initial";
    char solutionDatBaseFileName[MAX_FILE_NAME_LENGTH] = "solution";
    printf("Info: Parsing args\n");
    if (argc < 3 || argc > 5)
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
    char* endChr;
    int funcSelection = strtol(argv[2], &endChr, 10);
    if (endChr == argv[2] || funcSelection > 3 || funcSelection < 0)
    {
        PrintHelp();
        exit(1);
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

    if (argc > 4)
    {
        if(strlen(argv[4]) > MAX_FILE_NAME_LENGTH - 1)
        {
            printf("Error: Parameter file name too long\n");
            exit(1);
        }
        strncpy(solutionDatBaseFileName, argv[4], MAX_FILE_NAME_LENGTH);
    }

    GridParams params;
    if (ParseGridParameterFile(paramsFileName, &params))
    {
        printf("Error: Error in reading parameter file: %s\n", paramsFileName);
        exit(1);
    }
    printf("Info: Parameters:\n");
    PrintGridParameters(&params);
    printf("Info: Function selection: %d\n", funcSelection);

    const double dx = (params.m_XMax - params.m_XMin) / (params.m_NRow - 1);
    const double dy = (params.m_YMax - params.m_YMin) / (params.m_NCol - 1);

    double (*func)(double, double) = funcs[funcSelection];

    /* Produce both initial boundaries, and analytical solutions */
    double** gridInit = AllocateInitGrid(params.m_NRow, params.m_NCol);
    if (gridInit == NULL)
    {
        printf("Error: Error in allocating grid\n");
        exit(1);
    }
    double** gridSol = AllocateInitGrid(params.m_NRow, params.m_NCol);
    if (gridSol == NULL)
    {
        printf("Error: Error in allocating grid\n");
        /* Clean up */
        FreeGrid(gridInit);
        exit(1);
    }
    double x = params.m_XMin;
    double y = params.m_YMin;
    for (size_t i = 0; i < params.m_NRow; ++i)
    {
        y = params.m_YMin;
        for (size_t j = 0; j < params.m_NCol; ++j)
        {
            gridSol[i][j] = func(x, y);
            if (i > 0 && i < (params.m_NRow - 1) && j > 0 && j < (params.m_NCol - 1))
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
    snprintf(initialDatFileName, MAX_FILE_NAME_LENGTH, "%s.dat", initialDatBaseFileName);
    char solutionDatFileName[MAX_FILE_NAME_LENGTH];
    snprintf(solutionDatFileName, MAX_FILE_NAME_LENGTH, "%s.dat", solutionDatBaseFileName);

    WriteGrid(initialDatFileName, &params, gridInit);
    WriteGrid(solutionDatFileName, &params, gridSol);
    /* Clean up */
    FreeGrid(gridInit);
    FreeGrid(gridSol);
    printf("Info: Exiting\n");
}
