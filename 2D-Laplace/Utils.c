#include "Utils.h"

#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

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
    pprintf("Info: Parsing grid data from %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "r")) == NULL)
    {
        pprintf("Error: Cannot open %s for reading\n", fileName);
        return 1;
    }
    const double dx = (params->m_XMax - params->m_XMin) / (params->m_NRow - 1);
    const double dy = (params->m_YMax - params->m_YMin) / (params->m_NCol - 1);
    double x = params->m_XMin;
    double y = params->m_YMin;
    for (size_t i = 0; i < params->m_NRow; ++i)
    {
        y = params->m_YMin;
        for (size_t j = 0; j < params->m_NCol; ++j)
        {
            if(fscanf(fptr, "%*f %*f %lf\n", &grid1[i][j]) < 0)
            {
                pprintf("Error: Error in parsing %s\n", fileName);
                exit(1);
            }
            grid2[i][j] = grid1[i][j];
            //pprintf("READ: %f %f %f\n", x, y, grid1[i][j]);
            y += dy;
        }
        if(fscanf(fptr, "\n") < 0)
        {
            pprintf("Error: Error in parsing %s\n", fileName);
            exit(1);
        }
        x += dx;
    }
    if(fclose(fptr))
    {
        pprintf("Error: Cannot close %s\n", fileName);
        return 1;
    }

    return 0;
}

int WriteGrid(const char fileName[], const Params* params, double** grid)
{
    pprintf("Info: Writing grid data to %s\n", fileName);
    FILE* fptr;
    if ((fptr = fopen(fileName, "w")) == NULL)
    {
        pprintf("Error: Cannot open %s for writing\n", fileName);
        exit(1);
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
        pprintf("Error: Cannot close %s\n", fileName);
        exit(1);
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
    if (rank == 0)
    {
        va_list args;
        va_start(args, fmt);
        vprintf(fmt, args);
        va_end(args);
    }
    return 0;
}

