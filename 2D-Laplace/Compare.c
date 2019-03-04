/* *
 * Compare two .dat files
 */

#include "Utils.h"

#include <mpi.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define DEFAULT_TOLERANCE 0.001

void PrintHelp(void)
{
    printf("Usage: Compare <DatFile1> <DatFile2> [Tolerance = 1e-3] [Display Stats = 1]\n");
}

int main(int argc, char**argv)
{
    /* Parse args */
    double tolerance = DEFAULT_TOLERANCE;
    int displayStats = 1;
    printf("Info: Parsing args\n");
    if (argc < 3 || argc > 5)
    {
        PrintHelp();
        exit(1);
    }

    if (strlen(argv[1]) > MAX_FILE_NAME_LENGTH - 1)
    {
        printf("Error: File name too long\n");
        exit(1);
    }
    char datFileName1[MAX_FILE_NAME_LENGTH];
    strncpy(datFileName1, argv[1], MAX_FILE_NAME_LENGTH);

    if (strlen(argv[2]) > MAX_FILE_NAME_LENGTH - 1)
    {
        printf("Error: File name too long\n");
        exit(1);
    }
    char datFileName2[MAX_FILE_NAME_LENGTH];
    strncpy(datFileName2, argv[2], MAX_FILE_NAME_LENGTH);

    if (argc > 3)
    {
        char* endChr;
        tolerance = strtod(argv[3], &endChr);
        if (endChr == argv[3] || tolerance < 0)
        {
            PrintHelp();
            exit(1);
        }
    }

    if (argc > 4)
    {
        char* endChr;
        displayStats = strtol(argv[4], &endChr, 10);
        if (endChr == argv[4] || ! (displayStats == 0 || displayStats == 1))
        {
            PrintHelp();
            exit(1);
        }
    }

    GridParams gridParams1;
    if (ReadGridParams(datFileName1, &gridParams1))
    {
        exit(1);
    }
    GridParams gridParams2;
    if (ReadGridParams(datFileName2, &gridParams2))
    {
        exit(1);
    }
    if (!CompareGridParams(&gridParams1, &gridParams2))
    {
        printf("Error: Disagreeing parameters\n");
        exit(1);
    }

    int numMismatches = 0;

    double** grid1 = AllocateInitGrid(gridParams1.m_NRow, gridParams1.m_NCol);
    double** grid2 = AllocateInitGrid(gridParams2.m_NRow, gridParams2.m_NCol);
    if (ReadGrid(datFileName1, &gridParams1, grid1, NULL))
    {
        /* Clean up */
        FreeGrid(grid1);
        FreeGrid(grid2);
        exit(1);
    }
    if (ReadGrid(datFileName2, &gridParams2, grid2, NULL))
    {
        /* Clean up */
        FreeGrid(grid1);
        FreeGrid(grid2);
        exit(1);
    }

    for (size_t i = 0; i < gridParams1.m_NRow; ++i)
    {
        for (size_t j = 0; j < gridParams1.m_NCol; ++j)
        {
            if (fabs(grid1[i][j] - grid2[i][j]) > tolerance)
            {
                numMismatches++;
            }
        }
    }

    if (numMismatches > 0)
    {
        printf("Info: Mismatches: %d / %d (%f)\n", numMismatches, gridParams1.m_NRow * gridParams1.m_NCol,
                ((double) numMismatches )/ ((double)gridParams1.m_NRow * gridParams1.m_NCol));
        printf("Info: FAIL\n");
        return 1;
    }
    else
    {
        printf("Info: PASS\n");
        return 0;
    }

}
