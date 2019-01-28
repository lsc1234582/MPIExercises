#include "Utils.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* Analytical solutions
 * A series of functions that satisfy Laplace's equation
 */
double func0(double x, double y)
{
    return 1.0;
}

double func1(double x, double y)
{
    return 2 * x * y;
}

double func2(double x, double y)
{
    return x * x - y * y;
}

double func3(double x, double y)
{
    return pow(x, 4) - 6 * pow(x, 2) * pow(y, 2) + pow(y, 4);
}
/* Function map that holds the analytical solutions */
static double (*const funcs[])(double, double) =
{
    func0,
    func1,
    func2,
    func3
};

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
    char paramsFileName[128];
    if(strlen(argv[1]) > 128)
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
