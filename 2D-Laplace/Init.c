#include "Utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintHelp(void)
{
    printf("Usage: Init <ParamsFile> <FunctionSelection>\n\
            <FunctionSelect>: [0..3]\n");
}

int main(int argc, char**argv)
{
    /* Parse args */
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
    int funcSelection = strtol(argv[2], &endChr, 10);
    if (endChr == argv[2] || funcSelection > 3 || funcSelection < 0)
    {
        PrintHelp();
        exit(1);
    }

    Params params;
    if (ParseParameterFile(paramsFileName, &params))
    {
        printf("Error: Error in reading parameter file: %s\n", paramsFileName);
        exit(1);
    }
    printf("Info: Parameters:\n");
    PrintParameters(&params);
    printf("Info: Function selection: %d\n", funcSelection);

    const double dx = (params.m_XMax - params.m_XMin) / (params.m_NRow - 1);
    const double dy = (params.m_YMax - params.m_YMin) / (params.m_NCol - 1);

    double (*func)(double, double) = funcs[funcSelection];

    printf("Info: Writing initial bounadry values and analytical solutions to initial.dat and solution.dat\n");
    /* Produce both initial boundaries (writing to 'initial.dat'), and analytical solutions (writing to
     * 'solution.dat')) */
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

    if (WriteGrid("initial.dat", &params, gridInit))
    {
        exit(1);
    }
    if (WriteGrid("solution.dat", &params, gridSol))
    {
        exit(1);
    }
    /* Clean up */
    FreeGrid(params.m_NRow, gridInit);
    FreeGrid(params.m_NRow, gridSol);
    printf("Info: Exiting\n");
}
