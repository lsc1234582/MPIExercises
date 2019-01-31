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

int CreateParameterMPIStructDataType(MPI_Datatype* newType)
{
    const int numFields = 7;
    int blockLengths[] = {1, 1, 1, 1, 1, 1, 1};
    MPI_Datatype blockTypes[] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT, MPI_INT, MPI_DOUBLE};
    MPI_Aint blockDisplacements[numFields];
    blockDisplacements[0] = offsetof(Params, m_XMin);
    blockDisplacements[1] = offsetof(Params, m_XMax);
    blockDisplacements[2] = offsetof(Params, m_YMin);
    blockDisplacements[3] = offsetof(Params, m_YMax);
    blockDisplacements[4] = offsetof(Params, m_NRow);
    blockDisplacements[5] = offsetof(Params, m_NCol);
    blockDisplacements[6] = offsetof(Params, m_Tolerance);
    MPI_Type_create_struct(numFields, blockLengths, blockDisplacements, blockTypes, newType);
    MPI_Type_commit(newType);
    return 0;
}

int ParseParameterFile(const char fileName[], Params* params)
{
    FILE* fptr;
    if ((fptr = fopen(fileName, "r")) == NULL)
    {
        pprintf("Error: Cannot open parameter file %s for reading\n", fileName);
        return 1;
    }
    if(fscanf(fptr, "XMin: %lf\n", &params->m_XMin) < 0)
    {
        return 1;
    }
    if(fscanf(fptr, "XMax: %lf\n", &params->m_XMax) < 0)
    {
        return 1;
    }
    if(fscanf(fptr, "YMin: %lf\n", &params->m_YMin) < 0)
    {
        return 1;
    }
    if(fscanf(fptr, "YMax: %lf\n", &params->m_YMax) < 0)
    {
        return 1;
    }
    if(fscanf(fptr, "NRow: %d\n", &params->m_NRow) < 0)
    {
        return 1;
    }
    if(fscanf(fptr, "NCol: %d\n", &params->m_NCol) < 0)
    {
        return 1;
    }
    if(fscanf(fptr, "Tolerance: %lf\n", &params->m_Tolerance) < 0)
    {
        return 1;
    }
    return 0;
}

void PrintParameters(const Params* params)
{
    pprintf("XMin:\t%lf\nXMax:\t%lf\nYMin:\t%lf\nYMax:\t%lf\nNRow:\t%d\nNCol:\t%d\nTolerance:\t%lf\n",
            params->m_XMin, params->m_XMax, params->m_YMin, params->m_YMax,
            params->m_NRow, params->m_NCol, params->m_Tolerance);
}

double** AllocateInitGrid(const int nRow, const int nCol)
{
    double** grid = malloc(nRow * sizeof(double*));
    if (grid == NULL)
    {
        return NULL;
    }
    for (int i = 0; i < nRow; ++i)
    {
        if ((grid[i] = calloc(nCol, sizeof(double))) == NULL)
        {
            return NULL;
        }
    }
    return grid;
}

void FreeGrid(const int nRow, double** grid)
{
    for (int i = 0; i < nRow; ++i)
    {
        free(grid[i]);
    }
    free(grid);
}

int ReadGrid(const char fileName[], const Params* params, double** grid1, double** grid2)
{
    printf("Info: Parsing grid data from %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "r")) == NULL)
    {
        printf("Error: Cannot open %s for reading\n", fileName);
        fclose(fptr);
        return 1;
    }
    for (size_t i = 0; i < params->m_NRow; ++i)
    {
        for (size_t j = 0; j < params->m_NCol; ++j)
        {
            if(fscanf(fptr, "%*f %*f %lf\n", &grid1[i][j]) < 0)
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
        if(fscanf(fptr, "\n") < 0)
        {
            printf("Error: Error in parsing %s\n", fileName);
            fclose(fptr);
            return 1;
        }
    }
    fclose(fptr);
    return 0;
}

int ConcatenateGrid(const double** grid1, const double** grid2, const Params* param1, const Params* param2, const int axis, double** resultGrid, Params* resultParam)
{
    assert(axis == 0 || axis == 1);
    assert(param1->m_XMax <= param2->m_XMin);
    assert(param1->m_YMax <= param2->m_YMin);
    // Assert that dx and dy remain uniform before and after concatenating
    assert((param2->m_XMin - param1->m_XMax) == ((param1->m_XMax - param1->m_XMin) / param1->m_NCol) &&
           (param2->m_XMin - param1->m_XMax) == ((param2->m_XMax - param2->m_XMin) / param2->m_NCol));
    assert((param2->m_YMin - param1->m_YMax) == ((param1->m_YMax - param1->m_YMin) / param1->m_NRow) &&
           (param2->m_YMin - param1->m_YMax) == ((param2->m_YMax - param2->m_YMin) / param2->m_NRow));

    if (axis == 0)
    {
        assert(param1->m_NCol == param2->m_NCol);
        CopyGrid(grid1, param1, param1, resultGrid);
        CopyGrid(grid2, param2, param2, &resultGrid[param1->m_NRow]);
    }
    else if (axis == 1)
    {
        assert(param1->m_NRow == param2->m_NRow);
        for (int i = 0; i < param1->m_NRow + param2->m_NRow; ++i)
        {
            memcpy(&resultGrid[i][0], &grid1[i][0], param1->m_NCol * sizeof(double));
            memcpy(&resultGrid[i][param1->m_NCol], &grid2[i][0], param2->m_NCol * sizeof(double));
        }
    }
    resultParam->m_NCol = param1->m_NCol + param2->m_NCol;
    resultParam->m_NRow = param1->m_NRow + param2->m_NRow;
    resultParam->m_XMin = param1->m_XMin;
    resultParam->m_XMax = param2->m_XMax;
    resultParam->m_YMin = param1->m_YMin;
    resultParam->m_YMax = param2->m_YMax;
    return 0;
}

int CopyGrid(const double** srcGrid, const Params* srcParams, const Params* dstParams, double** dstGrid)
{
    assert(srcParams->m_NCol == dstParams->m_NCol && srcParams->m_NRow == dstParams->m_NRow);
    for (int i = 0; i < srcParams->m_NRow; ++i)
    {
        memcpy(dstGrid[i], srcGrid[i], srcParams->m_NCol * sizeof(double));
    }
    return 1;
}


int ConcatenateGridPatches(const double*** grids, const Params* params, const int numGridX, const int numGridY, double** resultGrid, Params* resultParam)
{
    int currPatchInx = 0;
    Params currFullPatchParam;
    double** currFullPatch;
    Params prevIncompleteFullPatchParam;
    double** prevIncompleteFullPatch;
    for (int i = 0; i < numGridX; ++i)
    {
        Params currHorPatchParam = {0.0, 0.0, 0.0, 0.0, params[currPatchInx].m_NRow, params[currPatchInx].m_NCol, 0.0};
        double** currHorPatch = AllocateInitGrid(currHorPatchParam.m_NRow, currHorPatchParam.m_NCol);
        CopyGrid(grids[currPatchInx], &currHorPatchParam, &currHorPatchParam, currHorPatch);
        Params prevIncompleteHorPatchParam = currHorPatchParam;
        double** prevIncompleteHorPatch = currHorPatch;
        for (int j = 0; j < numGridY - 1; ++j)
        {
            currPatchInx++;
            Params tempHorPatchParam = {0.0, 0.0, 0.0, 0.0, currHorPatchParam.m_NRow, currHorPatchParam.m_NCol + params[currPatchInx].m_NCol, 0.0};
            double** tempHorPatch = AllocateInitGrid(tempHorPatchParam.m_NRow, tempHorPatchParam.m_NCol);
            ConcatenateGrid((const double**)currHorPatch, grids[currPatchInx + 1], &currHorPatchParam, &params[currPatchInx + 1], 1, tempHorPatch, &tempHorPatchParam);
            currHorPatchParam = tempHorPatchParam;
            currHorPatch = tempHorPatch;
            FreeGrid(prevIncompleteHorPatchParam.m_NRow, prevIncompleteHorPatch);
            prevIncompleteHorPatchParam = currHorPatchParam;
            prevIncompleteHorPatch = currHorPatch;

        }
        if (i == 0)
        {
            currFullPatchParam = currHorPatchParam;
            currFullPatch = currHorPatch;
        }
        else
        {
            Params tempFullPatchParam = {0.0, 0.0, 0.0, 0.0, prevIncompleteFullPatchParam.m_NRow + currHorPatchParam.m_NRow, currHorPatchParam.m_NCol, 0.0};
            double** tempFullPatch = AllocateInitGrid(tempFullPatchParam.m_NRow, tempFullPatchParam.m_NCol);
            ConcatenateGrid((const double**)currFullPatch, (const double**)currHorPatch, &currFullPatchParam, &currHorPatchParam, 0, tempFullPatch, &tempFullPatchParam);
            currFullPatchParam = tempFullPatchParam;
            currFullPatch = tempFullPatch;
            FreeGrid(prevIncompleteFullPatchParam.m_NRow, prevIncompleteFullPatch);
            prevIncompleteFullPatchParam = currFullPatchParam;
            prevIncompleteFullPatch = currFullPatch;
        }
        prevIncompleteFullPatchParam = currFullPatchParam;
        prevIncompleteFullPatch = currFullPatch;
    }
    *resultParam = currFullPatchParam;
    resultGrid = currFullPatch;
    return 0;
}

/* A crude parser to read parameters from raw .dat file
 * Note that the format is more strict than the gnuplot .dat file in that each row of points must be followed by
 * exactly one empty line*/
int ReadGridParams(const char fileName[], Params* params)
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
            else if (sscanf(lineStr, "%[\n]", &tmp) == 1)
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
                printf("Error: Error in parsing %s\n, illegal characters %s\n", fileName, lineStr);
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
    PrintParameters(params);
    fclose(fptr);
    return 0;
}

int ReadGridHorPatch(const char fileName[], const Params* params, const GridHorPatch* horPatch, double** grid1, double** grid2)
{
    printf("Info: Parsing grid data from %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "r")) == NULL)
    {
        printf("Error: Cannot open %s for reading\n", fileName);
        fclose(fptr);
        return 1;
    }
    for (size_t i = 0; i < horPatch->m_NRow; ++i)
    {
        for (size_t j = 0; j < params->m_NCol; ++j)
        {
            if(fscanf(fptr, "%*f %*f %lf\n", &grid1[i][j]) < 0)
            {
                printf("Error: Error in parsing %s\n", fileName);
                fclose(fptr);
                return 1;
            }
            grid2[i][j] = grid1[i][j];
            //printf("READ: %f %f %f\n", x, y, grid1[i][j]);
        }
        if(fscanf(fptr, "\n") < 0)
        {
            printf("Error: Error in parsing %s\n", fileName);
            fclose(fptr);
            return 1;
        }
    }
    fclose(fptr);

    return 0;
}

int WriteGrid(const char fileName[], const Params* params, double** grid)
{
    printf("Info: Writing grid data to %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "w")) == NULL)
    {
        printf("Error: Cannot open %s for writing\n", fileName);
        return 1;
    }

    const double dx = (params->m_XMax - params->m_XMin) / (params->m_NRow - 1);
    const double dy = (params->m_YMax - params->m_YMin) / (params->m_NCol - 1);
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

int WriteGridHorPatch(const char fileName[], const Params* params, const GridHorPatch* horPatch, double** grid)
{
    printf("Info: Writing grid data to %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "w")) == NULL)
    {
        printf("Error: Cannot open %s for writing\n", fileName);
        return 1;
    }

    const double dx = (params->m_XMax - params->m_XMin) / (params->m_NRow - 1);
    const double dy = (params->m_YMax - params->m_YMin) / (params->m_NCol - 1);
    double x = horPatch->m_PatchX;
    double y = params->m_YMin;
    x = horPatch->m_PatchX;
    for (size_t i = 0; i < horPatch->m_NRow; ++i)
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

