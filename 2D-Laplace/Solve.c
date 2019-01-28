#include "Utils.h"
#include "Solvers.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintHelp(void)
{
    printf("Usage: Solve <ParamsFile>\n");
}

int main(int argc, char**argv)
{
    printf("Info: Parsing args\n");
    /* Parse args */
    if (argc != 2)
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
    Params params;

    if (ParseParameterFile(paramsFileName, &params))
    {
        printf("Error: Error in reading parameter file: %s\n", paramsFileName);
        exit(1);
    }
    printf("Info: Parameters:\n");
    PrintParameters(&params);

    SolveSerial(&params);

    printf("Info: Exiting\n");
}
