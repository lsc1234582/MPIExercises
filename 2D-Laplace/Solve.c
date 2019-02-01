#include "Utils.h"
#include "Solvers.h"

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintHelp(void)
{
    pprintf("Usage: Solve <ParamsFile>\n");
}

int main(int argc, char**argv)
{
    MPI_Init(&argc, &argv);
    pprintf("Info: Parsing args\n");
    /* Parse args */
    if (argc != 2)
    {
        PrintHelp();
        exit(1);
    }
    char paramsFileName[MAX_FILE_NAME_LENGTH];
    if(strlen(argv[1]) > MAX_FILE_NAME_LENGTH)
    {
        pprintf("Error: Parameter file name too long\n");
        exit(1);
    }
    strncpy(paramsFileName, argv[1], strlen(argv[1]));
    GridParams params;

    if (ParseGridParameterFile(paramsFileName, &params))
    {
        pprintf("Error: Error in reading parameter file: %s\n", paramsFileName);
        exit(1);
    }
    pprintf("Info: Parameters:\n");
    PrintGridParameters(&params);

    SolveSerial(&params);

    pprintf("Info: Exiting\n");
    MPI_Finalize();
}
