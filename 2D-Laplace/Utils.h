#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>

#define MASTER_RANK 0
#define MAX_FILE_NAME_LENGTH 128
#define EPSILON 0.001

/* Global Grid Parameters TODO: comments*/
typedef struct
{
    double m_XMin;
    double m_XMax;
    double m_YMin;
    double m_YMax;
    int m_NRow;
    int m_NCol;
    double m_Tolerance;
} GridParams;

/* Local Grid Patch Parameters TODO: comments*/
typedef struct
{
    int m_NRow; // Number of rows of the valid reading region.
    int m_NCol;
    int m_PatchI;
    int m_PatchJ;
    double m_PatchX; 
    double m_PatchY; 
    double m_Dx;
    double m_Dy;
    int m_LeftRank;
    int m_RightRank;
    int m_AboveRank;
    int m_BelowRank;
    int m_LeftMargin;  // Width of the 'halo' area
    int m_RightMargin;
    int m_AboveMargin;
    int m_BelowMargin;
    int m_NTotRow;   // Number of rows of the whole grid patch = m_NRow + m_LeftMargin + m_RightMargin
    int m_NTotCol;

    int m_LeftPadding;  // Width of the padding
    int m_RightPadding; 
    int m_AbovePadding;
    int m_BelowPadding;
} GridPatchParams;
/* Analytical solutions
 * A series of functions that satisfy Laplace's equation
 */
double func0(double x, double y);

double func1(double x, double y);

double func2(double x, double y);

double func3(double x, double y);

/* Function map that holds the analytical solutions */
extern double (*const funcs[])(double, double);

double GetDx(const GridParams* params);

double GetDy(const GridParams* params);

int CreateGridParameterMPIStructDataType(MPI_Datatype* newtype);

int CreateVertMarginMPIDataType(const GridPatchParams* patch, MPI_Datatype* newtype);

int ParseGridParameterFile(const char fileName[], GridParams* params);

void PrintGridParameters(const GridParams* params);

void PrintGrid(const GridParams* params, const double** grid);

void GetGridPatchParams(const GridParams* params, const int size, const int rank, const int nPatchInX, const int nPatchInY, GridPatchParams* patch);

double** AllocateInitGrid(const int nRow, const int nCol);

double** AllocateInitGridPatch(const GridPatchParams* patch);

void FreeGrid(double** grid);

int CopyGrid(const double** srcGrid, const GridParams* srcParams, const GridParams* dstParams, double** dstGrid);

int ReadGrid(const char fileName[], const GridParams* params, double** grid1, double** grid2);
int ReadGridPatch(const char fileName[], const GridPatchParams* patch, double** grid1, double** grid2);

int ReadGridParams(const char fileName[], GridParams* params);

int WriteGrid(const char fileName[], const GridParams* params, double** grid);
int WriteGridPatch(const char fileName[], const GridPatchParams* patch, double** grid);

int pprintf(const char* fmt, ...);

int ConcatenateGrid(const double** grid1, const double** grid2, const GridParams* param1, const GridParams* param2, const int axis, double** resultGrid, GridParams* resultParam);
int ConcatenateGrids(const double*** grids, const GridParams* params, const int numGridX, const int numGridY, double** resultGrid, GridParams* resultParam);


#endif
