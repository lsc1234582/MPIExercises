/* *
 * Utility functions
 */
#ifndef UTILS_H
#define UTILS_H

#include <mpi.h>

#define MASTER_RANK 0
#define MAX_FILE_NAME_LENGTH 128
#define EPSILON 0.001

/* Global Grid Parameters */
typedef struct
{
    double m_XMin;      // Minimum X value
    double m_XMax;      // Maximum X value
    double m_YMin;      // Minimum Y value
    double m_YMax;      // Maximum Y value
    int m_NRow;         // Number of rows (X axis)
    int m_NCol;         // Number of clumns (Y axis)
    double m_Tolerance; // Tolerance (convergence criteria)
} GridParams;

/* Local Grid Patch Parameters */
typedef struct
{
    int m_NRow;         // Number of rows of the Valid Region (X axis)
    int m_NCol;         // Number of rows of the Valid Region (Y axis)
    int m_NTotRow;      // Total number of rows = m_NRow + m_LeftMargin + m_RightMargin
    int m_NTotCol;      // Total number of columns = m_NRow + m_AboveMargin + m_BelowMargin

    int m_PatchI;       // The global row index of the first row (X axis)
    int m_PatchJ;       // The global column index of the first column (Y axis)
    double m_PatchX;    // The global X value of the first row (X axis)
    double m_PatchY;    // The global Y value of the first column (Y axis)
    double m_Dx;        // Delta x between rows
    double m_Dy;        // Delta y between columns
    int m_LeftRank;     // Rank of "left" process
    int m_RightRank;    // Rank of "right" process
    int m_AboveRank;    // Rank of "above" process
    int m_BelowRank;    // Rank of "below" process
    int m_LeftMargin;   // Width of the left side of the Halo Region
    int m_RightMargin;  // Width of the right side of the Halo Region
    int m_AboveMargin;  // Width of the top side of the Halo Region
    int m_BelowMargin;  // Width of the bottom side of the Halo Region

    int m_LeftPadding;  // Width of the left padding
    int m_RightPadding; // Width of the right padding
    int m_AbovePadding; // Width of the top padding
    int m_BelowPadding; // Width of the bottom padding
} GridPatchParams;

/* Analytical solutions
 * A series of arbitrary functions that satisfy Laplace's equation
 */
double func0(double x, double y);

double func1(double x, double y);

double func2(double x, double y);

double func3(double x, double y);

/* Function map that holds the analytical solutions */
extern double (*const funcs[])(double, double);

/* Calculate delta x */
double GetDx(const GridParams* params);

/* Calculate delta y */
double GetDy(const GridParams* params);

/* Create and commit an MPI data type for GridParams struct */
int CreateGridParameterMPIStructDataType(MPI_Datatype* newtype);

/* Create and commit an MPI data type for column margin elements */
int CreateColumnMarginElementMPIDataType(const GridPatchParams* patch, MPI_Datatype* newtype);

/* Initialise GridParams */
void InitialiseGridParams(GridParams* params);

/* Read GridParams from file */
int ParseGridParameterFile(const char fileName[], GridParams* params);

/* Allocate a local GridPatchParams */
void GetGridPatchParams(const GridParams* params, const int size, const int rank, const int nPatchInX, const int nPatchInY, GridPatchParams* patch);

/* Allocate and initialise a grid */
double** AllocateInitGrid(const int nRow, const int nCol);

/* Allocate and initialise a grid patch */
double** AllocateInitGridPatch(const GridPatchParams* patch);

/* Free a grid */
void FreeGrid(double** grid);

/* Copy a grid to another grid */
int CopyGrid(const double** srcGrid, const GridParams* srcParams, const GridParams* dstParams, double** dstGrid);

/* Read a grid from file */
int ReadGrid(const char fileName[], const GridParams* params, double** grid1, double** grid2);

/* Read a grid patch from file */
int ReadGridPatch(const char fileName[], const GridPatchParams* patch, double** grid1, double** grid2);

/* Compares two GridParams */
int CompareGridParams(const GridParams* params1, const GridParams* params2);

/* Read a GridParams from file */
int ReadGridParams(const char fileName[], GridParams* params);

/* Write a grid to file */
int WriteGrid(const char fileName[], const GridParams* params, double** grid);

/* Write a grid patch to file */
int WriteGridPatch(const char fileName[], const GridPatchParams* patch, double** grid);

/* Concate two grids in either x or y axis */
int ConcatenateGrid(const double** grid1, const double** grid2, const GridParams* param1, const GridParams* param2, const int axis, double** resultGrid, GridParams* resultParam);

/* Concatenate an array of grids */
int ConcatenateGrids(const double*** grids, const GridParams* params, const int numGridX, const int numGridY, double** resultGrid, GridParams* resultParam);

/* Print GridParams */
void PrintGridParameters(const GridParams* params);

/* Print Grid */
void PrintGrid(const GridParams* params, const double** grid);

/* Print command line arguments */
void PrintCLArgs(int argc, char** argv);

/* Parallel printf */
int pprintf(const char* fmt, ...);

#endif
