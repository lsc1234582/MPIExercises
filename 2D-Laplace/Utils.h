#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>

#define MASTER_RANK 0
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

/* Analytical solutions
 * A series of functions that satisfy Laplace's equation
 */
double func0(double x, double y);

double func1(double x, double y);

double func2(double x, double y);

double func3(double x, double y);

/* Function map that holds the analytical solutions */
extern double (*const funcs[])(double, double);

int CreateParameterMPIStructDataType(MPI_Datatype* newtype);

int ParseParameterFile(const char fileName[], Params* params);

void PrintParameters(const Params* params);

double** AllocateInitGrid(const int nRow, const int nCol);

void FreeGrid(const int nRow, double** grid);

int ReadGrid(const char fileName[], const Params* params, double** grid1, double** grid2);

int WriteGrid(const char fileName[], const Params* params, double** grid);

int pprintf(const char* fmt, ...);

#endif
