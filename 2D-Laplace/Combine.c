/* *
 * Combines multiple gnuplot .dat files into a single .dat file.
 */
#include "Utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintHelp(void)
{
    printf("Usage: Combine <NumPatchInX> <NumPatchInY> <DstDatFile> <SrcDatBaseFileName>\n");
}

int main(int argc, char** argv)
{
    /* Parse command-line artguments */
    if (argc != 5)
    {
        PrintHelp();
        exit(1);
    }
    char* endChr;
    int numPatchInX = strtol(argv[1], &endChr, 10);
    if (endChr == argv[1] || numPatchInX < 1)
    {
        PrintHelp();
        exit(1);
    }
    int numPatchInY = strtol(argv[2], &endChr, 10);
    if (endChr == argv[2] || numPatchInY < 1)
    {
        PrintHelp();
        exit(1);
    }
    int numPatches = numPatchInX * numPatchInY;
    char dstDatFileName[MAX_FILE_NAME_LENGTH];
    if(strlen(argv[3]) > MAX_FILE_NAME_LENGTH - 1)
    {
        printf("Error: .dat file name too long\n");
        exit(1);
    }
    strncpy(dstDatFileName, argv[3], MAX_FILE_NAME_LENGTH);
    char srcDatFileBaseName[MAX_FILE_NAME_LENGTH];
    if(strlen(argv[4]) > MAX_FILE_NAME_LENGTH - 1)
    {
        printf("Error: .dat file name too long\n");
        exit(1);
    }
    strncpy(srcDatFileBaseName, argv[4], MAX_FILE_NAME_LENGTH);

    /* Read source grid parameters and source grid data */
    GridParams params[numPatches];
    double** grids[numPatches];
    for (int i = 0; i < numPatches; ++i)
    {
        char srcDatFileName[MAX_FILE_NAME_LENGTH];
        snprintf(srcDatFileName, MAX_FILE_NAME_LENGTH, "%s.MPI_%d.dat", srcDatFileBaseName, i);
        ReadGridParams(srcDatFileName, &params[i]);
        grids[i] = AllocateInitGrid(params[i].m_NRow, params[i].m_NCol);
        if (ReadGrid(srcDatFileName, &params[i], grids[i], NULL))
        {
            /* Clean up */
            for (int j = 0; j <= i; ++j)
            {
                FreeGrid(grids[j]);
            }
            exit(1);
        }
    }

    int totalNumRow = 0;
    int totalNumCol = 0;
    for (int i = 0; i < numPatchInX; ++i)
    {
        totalNumRow += params[i * numPatchInY].m_NRow;
    }
    for (int j = 0; j < numPatchInY; ++j)
    {
        totalNumCol += params[j].m_NCol;
    }

    /* Perform concatenation */
    GridParams dstParams;
    double** dstGrid = AllocateInitGrid(totalNumRow, totalNumCol);
    ConcatenateGrids((const double***)grids, params, numPatchInX, numPatchInY, dstGrid, &dstParams);

    /* Write destination grid data */
    WriteGrid(dstDatFileName, &dstParams, dstGrid);

    /* Clean up */
    for (int i = 0; i < numPatches; ++i)
    {
        FreeGrid(grids[i]);
    }
    FreeGrid(dstGrid);
    printf("Info: Exiting\n");
    return 0;
}
