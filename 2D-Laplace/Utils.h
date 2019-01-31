#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>

#define MASTER_RANK 0
#define MAX_FILE_NAME_LENGTH 128
#define EPSILON 0.00001
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


typedef struct
{
    int m_NRow;
    int m_PatchI;
    double m_PatchX; 
} GridHorPatch;
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

int CopyGrid(const double** srcGrid, const Params* srcParams, const Params* dstParams, double** dstGrid);

int ReadGrid(const char fileName[], const Params* params, double** grid1, double** grid2);
int ReadGridHorPatch(const char fileName[], const Params* params, const GridHorPatch* horPatch, double** grid1, double** grid2);

int ReadGridParams(const char fileName[], Params* params);

int WriteGrid(const char fileName[], const Params* params, double** grid);
int WriteGridHorPatch(const char fileName[], const Params* params, const GridHorPatch* horPatch, double** grid);

int pprintf(const char* fmt, ...);

int ConcatenateGrid(const double** grid1, const double** grid2, const Params* param1, const Params* param2, const int axis, double** resultGrid, Params* resultParam);

int ConcatenateGridPatches(const double*** grids, const Params* params, const int numGridX, const int numGridY, double** resultGrid, Params* resultParam);
#endif
