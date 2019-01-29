#ifndef UTILS_H
#define UTILS_H

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

int ParseParameterFile(const char fileName[], Params* params);

void PrintParameters(const Params* params);

double** AllocateInitGrid(const int nRow, const int nCol);

void FreeGrid(const int nRow, double** grid);

int ReadGrid(const char fileName[], const Params* params, double** grid1, double** grid2);

int WriteGrid(const char fileName[], const Params* params, double** grid);

int pprintf(const char* fmt, ...);

#endif
