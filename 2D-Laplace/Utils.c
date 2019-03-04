#include "Utils.h"

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
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
double (*const funcs[])(double, double) =
{
    func0,
    func1,
    func2,
    func3
};

void InitialiseGridParams(GridParams* params)
{
    params->m_XMin = 0.0;
    params->m_XMax = 0.0;
    params->m_YMin = 0.0;
    params->m_YMax = 0.0;
    params->m_NRow = 0;
    params->m_NCol = 0;
    params->m_Tolerance = -1.0;
}

double GetDx(const GridParams* params)
{
    assert(params->m_XMax > params->m_XMin);
    assert(params->m_NRow > 1);
    return (params->m_XMax - params->m_XMin) / (params->m_NRow - 1);
}

double GetDy(const GridParams* params)
{
    assert(params->m_YMax > params->m_YMin);
    assert(params->m_NCol > 1);
    return (params->m_YMax - params->m_YMin) / (params->m_NCol - 1);
}

int CreateGridParameterMPIStructDataType(MPI_Datatype* newType)
{
    const int numFields = 7;
    int blockLengths[] = {1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype blockTypes[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_DOUBLE};
    MPI_Aint blockDisplacements[numFields];
    blockDisplacements[0] = offsetof(GridParams, m_XMin);
    blockDisplacements[1] = offsetof(GridParams, m_XMax);
    blockDisplacements[2] = offsetof(GridParams, m_YMin);
    blockDisplacements[3] = offsetof(GridParams, m_YMax);
    blockDisplacements[4] = offsetof(GridParams, m_NRow);
    blockDisplacements[5] = offsetof(GridParams, m_NCol);
    blockDisplacements[6] = offsetof(GridParams, m_Tolerance);
    MPI_Type_create_struct(numFields, blockLengths, blockDisplacements, blockTypes, newType);
    MPI_Type_commit(newType);
    return 0;
}

int CreateColumnMarginElementMPIDataType(const GridPatchParams* patch, MPI_Datatype* newType)
{
    MPI_Type_create_resized(MPI_DOUBLE, 0, (MPI_Aint)(patch->m_NTotCol * sizeof(MPI_DOUBLE)), newType);
    MPI_Type_commit(newType);
    return 1;
}

int ParseGridParameterFile(const char fileName[], GridParams* params)
{
    FILE* fptr;
    if ((fptr = fopen(fileName, "r")) == NULL)
    {
        pprintf("Error: Cannot open parameter file %s for reading\n", fileName);
        return 1;
    }
    if(fscanf(fptr, "XMin: %lf XMax: %lf YMin: %lf YMax: %lf NRow: %d NCol: %d Tolerance: %lf",
                &params->m_XMin,
                &params->m_XMax,
                &params->m_YMin,
                &params->m_YMax,
                &params->m_NRow,
                &params->m_NCol,
                &params->m_Tolerance) < 7)
    {
        printf("Error: Error in parsing parameter file %s\n", fileName);
        fclose(fptr);
        return 1;
    }
    fclose(fptr);
    return 0;
}

void PrintGridParameters(const GridParams* params)
{
    pprintf("XMin:\t%lf\nXMax:\t%lf\nYMin:\t%lf\nYMax:\t%lf\nNRow:\t%d\nNCol:\t%d\nTolerance:\t%lf\n",
            params->m_XMin, params->m_XMax, params->m_YMin, params->m_YMax,
            params->m_NRow, params->m_NCol, params->m_Tolerance);
}

/* Allocate a continuous block of memory for a 2D array (grid) */
double** AllocateInitGrid(const int nRow, const int nCol)
{
    double** grid = (double**) malloc(nRow * sizeof(double*));
    if (grid == NULL)
    {
        return NULL;
    }
    grid[0] = (double*) malloc(nRow * nCol * sizeof(double));
    if (grid[0] == NULL)
    {
        return NULL;
    }
    for (int i = 1; i < nRow; ++i)
    {
        grid[i] = grid[i - 1] + nCol;
    }
    return grid;
}

double** AllocateInitGridPatch(const GridPatchParams* patch)
{
    return AllocateInitGrid(patch->m_NTotRow, patch->m_NTotCol);
}

void FreeGrid(double** grid)
{
    free(grid[0]);
    free(grid);
}

int ConcatenateGrid(const double** grid1, const double** grid2, const GridParams* param1, const GridParams* param2, const int axis, double** resultGrid, GridParams* resultParam)
{
    assert(axis == 0 || axis == 1);

    if (axis == 0)
    {
        assert(param1->m_NCol == param2->m_NCol);
        assert(param1->m_XMax <= param2->m_XMin);
        // Assert that dx is uniform before concatenating
        double dxDiff = fabs(GetDx(param1) - GetDx(param2));
        if (dxDiff > EPSILON)
        {
            printf("Warning: significant dx difference: %lf\n", dxDiff);
        }
        CopyGrid(grid1, param1, param1, resultGrid);
        CopyGrid(grid2, param2, param2, &resultGrid[param1->m_NRow]);
        resultParam->m_NRow = param1->m_NRow + param2->m_NRow;
        resultParam->m_XMin = param1->m_XMin;
        resultParam->m_XMax = param2->m_XMax;
        // Assert that dx remains uniform after concatenating
        dxDiff = fabs(GetDx(param1) - GetDx(resultParam));
        if (dxDiff > EPSILON)
        {
            printf("Warning: significant dx difference: %lf\n", dxDiff);
        }
    }
    else if (axis == 1)
    {
        assert(param1->m_NRow == param2->m_NRow);
        //printf("param1->m_YMax: %f, param2->m_YMin: %f\n", param1->m_YMax, param2->m_YMin);
        assert(param1->m_YMax <= param2->m_YMin);
        // Assert that dy is uniform before concatenating
        double dyDiff = fabs(GetDy(param1) - GetDy(param2));
        if (dyDiff > EPSILON)
        {
            printf("Warning: significant dy difference: %lf\n", dyDiff);
        }
        for (int i = 0; i < param1->m_NRow; ++i)
        {
            memcpy(&resultGrid[i][0], &grid1[i][0], param1->m_NCol * sizeof(double));
            memcpy(&resultGrid[i][param1->m_NCol], &grid2[i][0], param2->m_NCol * sizeof(double));
        }
        resultParam->m_NCol = param1->m_NCol + param2->m_NCol;
        resultParam->m_YMin = param1->m_YMin;
        resultParam->m_YMax = param2->m_YMax;
        // Assert that dy remains uniform after concatenating
        dyDiff = fabs(GetDy(param1) - GetDy(resultParam));
        if (dyDiff > EPSILON)
        {
            printf("Warning: significant dy difference: %lf\n", dyDiff);
        }
    }
    else
    {
        // Unsupported;
        assert(0);
    }
    return 0;
}

int CopyGrid(const double** srcGrid, const GridParams* srcParams, const GridParams* dstParams, double** dstGrid)
{
    assert(srcParams->m_NCol == dstParams->m_NCol && srcParams->m_NRow == dstParams->m_NRow);
    for (int i = 0; i < srcParams->m_NRow; ++i)
    {
        memcpy(dstGrid[i], srcGrid[i], srcParams->m_NCol * sizeof(double));
    }
    return 0;
}

void PrintGrid(const GridParams* params, const double** grid)
{
    for (int i = 0; i < params->m_NRow; ++i)
    {
        for (int j = 0; j < params->m_NCol; ++j)
        {
            printf("%lf ", grid[i][j]);
        }
        printf("\n");
    }
}

int ConcatenateGrids(const double*** grids, const GridParams* params, const int numGridX, const int numGridY, double** resultGrid, GridParams* resultParam)
{
    int currPatchInx = 0;
    GridParams currFullPatchParam;
    double** currFullPatch;
    GridParams prevIncompleteFullPatchParam;
    double** prevIncompleteFullPatch;
    for (int i = 0; i < numGridX; ++i)
    {
        GridParams currHorPatchParam = params[currPatchInx];
        double** currHorPatch = AllocateInitGrid(currHorPatchParam.m_NRow, currHorPatchParam.m_NCol);
        CopyGrid(grids[currPatchInx], &currHorPatchParam, &currHorPatchParam, currHorPatch);
        GridParams prevIncompleteHorPatchParam = currHorPatchParam;
        double** prevIncompleteHorPatch = currHorPatch;
        for (int j = 0; j < numGridY; ++j)
        {
            if (j > 0)
            {
                GridParams tempHorPatchParam = currHorPatchParam;
                tempHorPatchParam.m_NCol += params[currPatchInx].m_NCol;
                double** tempHorPatch = AllocateInitGrid(tempHorPatchParam.m_NRow, tempHorPatchParam.m_NCol);
                ConcatenateGrid((const double**)currHorPatch, grids[currPatchInx], &currHorPatchParam, &params[currPatchInx], 1, tempHorPatch, &tempHorPatchParam);
                currHorPatchParam = tempHorPatchParam;
                currHorPatch = tempHorPatch;
                FreeGrid(prevIncompleteHorPatch);
                prevIncompleteHorPatchParam = currHorPatchParam;
                prevIncompleteHorPatch = currHorPatch;
            }
            currPatchInx++;
        }
        if (i == 0)
        {
            currFullPatchParam = currHorPatchParam;
            currFullPatch = currHorPatch;
        }
        else
        {
            GridParams tempFullPatchParam = currFullPatchParam;
            tempFullPatchParam.m_NRow += currHorPatchParam.m_NRow;
            double** tempFullPatch = AllocateInitGrid(tempFullPatchParam.m_NRow, tempFullPatchParam.m_NCol);
            ConcatenateGrid((const double**)currFullPatch, (const double**)currHorPatch, &currFullPatchParam, &currHorPatchParam, 0, tempFullPatch, &tempFullPatchParam);
            currFullPatchParam = tempFullPatchParam;
            currFullPatch = tempFullPatch;
            FreeGrid(prevIncompleteFullPatch);
            prevIncompleteFullPatchParam = currFullPatchParam;
            prevIncompleteFullPatch = currFullPatch;
        }
        prevIncompleteFullPatchParam = currFullPatchParam;
        prevIncompleteFullPatch = currFullPatch;
    }
    *resultParam = currFullPatchParam;
    // Perform a shallow copy
    for (int i = 0; i < resultParam->m_NRow; ++i)
    {
        resultGrid[i] = currFullPatch[i];
    }
    return 0;
}

int CompareGridParams(const GridParams* params1, const GridParams* params2)
{
    int result = 1;
    result &= fabs(params1->m_XMin - params2->m_XMin) < EPSILON;
    result &= fabs(params1->m_XMax - params2->m_XMax) < EPSILON;
    result &= fabs(params1->m_YMin - params2->m_YMin) < EPSILON;
    result &= fabs(params1->m_YMax - params2->m_YMax) < EPSILON;
    result &= params1->m_NRow == params2->m_NRow;
    result &= params1->m_NCol == params2->m_NCol;
    result &= fabs(params1->m_Tolerance - params2->m_Tolerance) < EPSILON;
    return result;
}

/* A crude parser to read parameters from raw .dat file
 * Note that the format is more strict than the gnuplot .dat file in that each row of points must be followed by
 * exactly one empty line*/
int ReadGridParams(const char fileName[], GridParams* params)
{
    printf("Info: Parsing grid parameters from %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "r")) == NULL)
    {
        printf("Error: Cannot open %s for reading\n", fileName);
        fclose(fptr);
        return 1;
    }

    char lineStr[256];
    char tmp;
    int i = 0;
    int j = 0;
    double lx = 0.0;
    double ly = 0.0;
    double x = 0.0;
    double y = 0.0;
    double ldx, ldy, dx, dy;
    params->m_NRow = -1;
    params->m_NCol = -1;
    while (1)
    {
        while (fgets(lineStr, 256, fptr) != NULL)
        {
            if (sscanf(lineStr, "%lf %lf %*f", &x, &y) == 2)
            {
                //printf("%d, %d, %f, %f\n", i, j, x, y);
                if (i == 0 && j == 0)
                {
                    params->m_XMin = x;
                    params->m_YMin = y;
                }
                if (i > 0)
                {
                    dx = x - lx;
                    if (dx < 0)
                    {
                        printf("Error: not standard coordiate system (expecting starting with X min)\n");
                        fclose(fptr);
                        return 1;
                    }
                    if (i > 1 && fabs(ldx - dx) > EPSILON)
                    {
                        printf("Warning: non-uniform dx: %lf %lf\n", ldx, dx);
                    }
                    ldx = dx;
                }
                if (j > 0)
                {
                    dy = y - ly;
                    if (dy < 0)
                    {
                        printf("Error: not standard coordiate system (expecting starting with Y min)\n");
                        fclose(fptr);
                        return 1;
                    }
                    if (j > 1 && fabs(ldy - dy) > EPSILON)
                    {
                        printf("Warning: non-uniform dy: %lf %lf\n", ldy, dy);
                    }
                    ldy = dy;
                }
                ly = y;
                j++;
            }
            else if (sscanf(lineStr, "%1c", &tmp) == 1 && tmp == '\n')
            {
                //printf("New line\n");
                if (params->m_NCol == -1)
                {
                    params->m_NCol = j;
                }
                else if (params->m_NCol != j)
                {
                    printf("Error: Error in parsing %s, nonmatching dimensions: %d %d\n", fileName, params->m_NCol, j);
                    fclose(fptr);
                    return 1;
                }
                j = 0;
                i++;
                lx = x;
            }
            else
            {
                printf("Error: Error in parsing %s, illegal characters %s\n", fileName, lineStr);
                fclose(fptr);
                return 1;
            }
        }
        if (feof(fptr))
        {
            break;
        }
        else
        {
            printf("Error: Error in parsing %s\n", fileName);
            fclose(fptr);
            return 1;
        }
    }
    params->m_XMax = x;
    params->m_YMax = y;
    params->m_NRow = i;
#ifdef DEBUG
    PrintGridParameters(params);
#endif
    fclose(fptr);
    return 0;
}

int ReadGrid(const char fileName[], const GridParams* params, double** grid1, double** grid2)
{
    printf("Info: Parsing grid data from %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "r")) == NULL)
    {
        printf("Error: Cannot open %s for reading\n", fileName);
        return 1;
    }
    for (size_t i = 0; i < params->m_NRow; ++i)
    {
        for (size_t j = 0; j < params->m_NCol; ++j)
        {
            if(fscanf(fptr, "%*f %*f %lf", &grid1[i][j]) < 0)
            {
                printf("Error: Error in parsing %s\n", fileName);
                fclose(fptr);
                return 1;
            }
            if (grid2 != NULL)
            {
                grid2[i][j] = grid1[i][j];
            }
            //printf("READ: %f %f %f\n", x, y, grid1[i][j]);
        }
    }
    fclose(fptr);
    return 0;
}

int ReadGridPatch(const char fileName[], const GridPatchParams* patch, double** grid1, double** grid2)
{
    printf("Info: Parsing grid data from %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "r")) == NULL)
    {
        printf("Error: Cannot open %s for reading\n", fileName);
        return 1;
    }
    for (size_t i = patch->m_AboveMargin; i < patch->m_NTotRow - patch->m_BelowMargin; ++i)
    {
        for (size_t j = patch->m_LeftMargin; j < patch->m_NTotCol - patch->m_RightMargin; ++j)
        {
            if(fscanf(fptr, "%*f %*f %lf", &grid1[i][j]) < 0)
            {
                printf("Error: Error in parsing %s\n", fileName);
                fclose(fptr);
                return 1;
            }
            if (grid2 != NULL)
            {
                grid2[i][j] = grid1[i][j];
            }
            //printf("READ: %f %f %f\n", x, y, grid1[i][j]);
        }
    }
    fclose(fptr);
    return 0;
}

int WriteGrid(const char fileName[], const GridParams* params, double** grid)
{
    printf("Info: Writing grid data to %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "w")) == NULL)
    {
        printf("Error: Cannot open %s for writing\n", fileName);
        return 1;
    }

    const double dx = GetDx(params);
    const double dy = GetDy(params);
    double x = params->m_XMin;
    double y = params->m_YMin;
    x = params->m_XMin;
    for (size_t i = 0; i < params->m_NRow; ++i)
    {
        y = params->m_YMin;
        for (size_t j = 0; j < params->m_NCol; ++j)
        {
            fprintf(fptr, "%f %f %f\n", x, y, grid[i][j]);
            y += dy;
        }
        fprintf(fptr, "\n");
        x += dx;
    }
    if(fclose(fptr))
    {
        printf("Error: Cannot close %s\n", fileName);
        return 1;
    }
    return 0;
}

int WriteGridPatch(const char fileName[], const GridPatchParams* patch, double** grid)
{
    printf("Info: Writing grid data to %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "w")) == NULL)
    {
        printf("Error: Cannot open %s for writing\n", fileName);
        return 1;
    }

    double x = patch->m_PatchX;
    double y = patch->m_PatchY;
    x = patch->m_PatchX;
    for (size_t i = patch->m_AboveMargin; i < patch->m_NTotRow - patch->m_BelowMargin; ++i)
    {
        y = patch->m_PatchY;
        for (size_t j = patch->m_LeftMargin; j < patch->m_NTotCol - patch->m_RightMargin; ++j)
        {
            fprintf(fptr, "%f %f %f\n", x, y, grid[i][j]);
            y += patch->m_Dy;
        }
        fprintf(fptr, "\n");
        x += patch->m_Dx;
    }
    if(fclose(fptr))
    {
        printf("Error: Cannot close %s\n", fileName);
        return 1;
    }
    return 0;
}

int pprintf(const char* fmt, ...)
{
    int initialised;
    int rank;
    if (MPI_Initialized(&initialised) || !initialised)
    {
        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
        return 1;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == MASTER_RANK)
    {
        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
    }
    return 0;
}

static int CoordToIndColMajor(const int rowInd, const int colInd, const int nRow, const int nCol)
{
    assert(nRow > 0 && nCol > 0);
    assert(rowInd >= 0 && rowInd < nRow);
    assert(colInd >= 0 && colInd < nCol);
    return rowInd * nCol + colInd;
}

static void IndToCoordColMajor(const int ind, const int nRow, const int nCol, int* rowInd, int* colInd)
{
    assert(nRow > 0 && nCol > 0);
    assert(ind >= 0 && ind < nRow * nCol);
    *colInd = ind % nCol;
    *rowInd = floor(ind / nCol);
}

void GetGridPatchParams(const GridParams* params, const int size, const int rank, const int nPatchInX, const int nPatchInY, GridPatchParams* patch)
{
    assert(size == nPatchInX * nPatchInY);
    const int fullRowSize = (int) ceil((double)params->m_NRow / (double)nPatchInX);
    const int partialRowSize = params->m_NRow % fullRowSize;
    const int fullColSize = (int) ceil((double)params->m_NCol / (double)nPatchInY);
    const int partialColSize = params->m_NCol % fullColSize;
    const double dx = GetDx(params);
    const double dy = GetDy(params);

    int patchI, patchJ;
    IndToCoordColMajor(rank, nPatchInX, nPatchInY, &patchI, &patchJ);

    patch->m_NRow = (partialRowSize == 0 || patchI < nPatchInX - 1) ? fullRowSize : partialRowSize;
    patch->m_NCol = (partialColSize == 0 || patchJ < nPatchInY - 1) ? fullColSize : partialColSize;
    patch->m_PatchI = patchI * fullRowSize;
    patch->m_PatchJ = patchJ * fullColSize;
    patch->m_PatchX = params->m_XMin + patch->m_PatchI * dx;
    patch->m_PatchY = params->m_YMin + patch->m_PatchJ * dy;
    patch->m_Dx = dx;
    patch->m_Dy = dy;
    patch->m_AboveRank = patchI > 0 ? CoordToIndColMajor(patchI - 1, patchJ, nPatchInX, nPatchInY) : MPI_PROC_NULL;
    patch->m_BelowRank = patchI < nPatchInX - 1 ? CoordToIndColMajor(patchI + 1, patchJ, nPatchInX, nPatchInY) : MPI_PROC_NULL;
    patch->m_LeftRank = patchJ > 0 ? CoordToIndColMajor(patchI, patchJ - 1, nPatchInX, nPatchInY) : MPI_PROC_NULL;
    patch->m_RightRank = patchJ < nPatchInY - 1 ? CoordToIndColMajor(patchI, patchJ + 1, nPatchInX, nPatchInY) : MPI_PROC_NULL;
    // Do not have any 'halo'/margin in x/y if there's no division happening in that axis (nPatchIn... == 1) or if the
    // patch is located at the edge
    patch->m_AboveMargin = (nPatchInX > 1 && patchI > 0) ? 1 : 0;
    patch->m_BelowMargin = (nPatchInX > 1 && patchI < nPatchInX - 1) ? 1 : 0;
    patch->m_LeftMargin = (nPatchInY > 1 && patchJ > 0) ? 1 : 0;
    patch->m_RightMargin = (nPatchInY > 1 && patchJ < nPatchInY - 1) ? 1 : 0;
    patch->m_NTotRow = patch->m_AboveMargin + patch->m_NRow + patch->m_BelowMargin;
    patch->m_NTotCol = patch->m_LeftMargin + patch->m_NCol + patch->m_RightMargin;
    patch->m_AbovePadding = 1;
    patch->m_BelowPadding = 1;
    patch->m_LeftPadding = 1;
    patch->m_RightPadding = 1;
}

void PrintCLArgs(int argc, char** argv)
{
    for (int i = 0; i < argc; ++i)
    {
        pprintf("%s ", argv[i]);
    }
    pprintf("\n");
}
