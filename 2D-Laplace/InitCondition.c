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
    printf("Usage: InitCondition <ParamsFile> <FunctionSelection>\n\
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

    FILE* fptrInit;
    FILE* fptrSol;
    if ((fptrInit = fopen("initial.dat", "w")) == NULL)
    {
        printf("Error: Cannot open initial.dat for writing\n");
        exit(1);
    }

    if ((fptrSol = fopen("solution.dat", "w")) == NULL)
    {
        printf("Error: Cannot open solution.dat for writing\n");
        exit(1);
    }
    double (*func)(double, double) = funcs[funcSelection];

    printf("Info: Writing initial bounadry values and analytical solutions to initial.dat and solutions.dat\n");
    /* Produce both initial boundaries (writing to 'initial.dat'), and analytical solutions (writing to
     * 'solution.dat')) */
    double x = params.m_XMin;
    double y = params.m_YMin;
    for (size_t i = 0; i < params.m_NRow; ++i)
    {
        y = params.m_YMin;
        for (size_t j = 0; j < params.m_NCol; ++j)
        {
            if (i > 0 && i < (params.m_NRow - 1) && j > 0 && j < (params.m_NCol - 1))
            {
                fprintf(fptrInit, "%f %f %f\n", x, y, 0.0);
            }
            else
            {
                fprintf(fptrInit, "%f %f %f\n", x, y, func(x, y));
            }
            fprintf(fptrSol, "%f %f %f\n", x, y, func(x, y));
            y += dy;
        }
        fprintf(fptrInit, "\n");
        fprintf(fptrSol, "\n");
        x += dx;
    }

    if(fclose(fptrInit))
    {
        printf("Error: Cannot close initial.dat\n");
        exit(1);
    }
    if(fclose(fptrSol))
    {
        printf("Error: Cannot close solution.dat\n");
        exit(1);
    }
    printf("Info: Exiting\n");
}
