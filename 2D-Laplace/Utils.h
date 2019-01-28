#ifndef UTILS_H
#define UTILS_H
#include <stdio.h>

/* Parameters TODO: comments*/
typedef struct
{
    double m_XMin;
    double m_XMax;
    double m_YMin;
    double m_YMax;
    int m_NRow;
    int m_NCol;
    double m_Tolerance;
} Params;

int ParseParameterFile(const char fileName[], Params* params)
{
    FILE* fptr;
    if ((fptr = fopen(fileName, "r")) == NULL)
    {
        printf("Error: Cannot open parameter file %s for reading\n", fileName);
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
    printf("XMin:\t%lf\nXMax:\t%lf\nYMin:\t%lf\nYMax:\t%lf\nNRow:\t%d\nNCol:\t%d\nTolerance:\t%lf\n",
            params->m_XMin, params->m_XMax, params->m_YMin, params->m_YMax,
            params->m_NRow, params->m_NCol, params->m_Tolerance);
}

#endif
