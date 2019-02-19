#include "Utils.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void PrintHelp(void)
{
    printf("Usage: Combine <NumPatchInX> <NumPatchInY> <DstDatFile> <SrcDatFile0>...\n");
}

int main(int argc, char** argv)
{
    if (argc < 5)
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
    if (numPatches != argc - 4)
    {
        printf("Error: Number of patches mismatches number of provided .dat files\n");
        exit(1);
    }
    char dstDatFileName[MAX_FILE_NAME_LENGTH];
    if(strlen(argv[3]) > MAX_FILE_NAME_LENGTH)
    {
        printf("Error: .dat file name too long\n");
        exit(1);
    }
    strncpy(dstDatFileName, argv[3], strlen(argv[3]));
    char srcDatFileNames[numPatches][MAX_FILE_NAME_LENGTH];
    for (int i = 4; i < numPatches + 4; ++i)
    {
        if(strlen(argv[i]) > MAX_FILE_NAME_LENGTH)
        {
            printf("Error: .dat file name too long\n");
            exit(1);
        }
        strncpy(srcDatFileNames[i - 4], argv[i], MAX_FILE_NAME_LENGTH);
        //printf("FILE: %s\n", srcDatFileNames[i - 4]);
    }

    /* Read source grid parameters and source grid data */
    GridParams params[numPatches];
    double** grids[numPatches];
    for (int i = 0; i < numPatches; ++i)
    {
        ReadGridParams(srcDatFileNames[i], &params[i]);
        grids[i] = AllocateInitGrid(params[i].m_NRow, params[i].m_NCol);
        ReadGrid(srcDatFileNames[i], &params[i], grids[i], NULL);
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
        free(grids[i]);
    }
    free(dstGrid);
    return 0;
}
