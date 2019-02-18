#include "Solvers.h"

#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

int SolveSerial(const GridParams* params)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {
        /* Initialise grid by reading "initial.dat"
         * Set boundary values to those set by the analitycal solution (initial.dat), and all the interior values to zeros.
         */
        pprintf("Info: Solving 2d Laplace with serial Jacobi iteration method\n");
        double** grid1 = AllocateInitGrid(params->m_NRow, params->m_NCol);
        if (grid1 == NULL)
        {
            pprintf("Error: Error in allocating grid\n");
            exit(1);
        }
        double** grid2 = AllocateInitGrid(params->m_NRow, params->m_NCol);
        if (grid2 == NULL)
        {
            pprintf("Error: Error in allocating grid\n");
            exit(1);
        }


        /* Parse initial.dat */
        if (ReadGrid("initial.dat", params, grid1, grid2))
        {
            exit(1);
        }

        /* Solve boundary value problem with Jacobi iteration method */
        pprintf("Info: Solving...\n");
        const double dx = GetDx(params);
        const double dy = GetDy(params);
        double x = params->m_XMin;
        double y = params->m_YMin;
        double maxDiff;
        double** tempGrid = NULL;
        do
        {
            maxDiff = 0.0;
            x = params->m_XMin;
            for (size_t i = 1; i < params->m_NRow - 1; ++i)
            {
                y = params->m_YMin;
                for (size_t j = 1; j < params->m_NCol - 1; ++j)
                {
                    double term1 = (grid1[i-1][j] + grid1[i+1][j]) / (dx * dx) + (grid1[i][j-1] + grid1[i][j+1]) / (dy * dy); 
                    double term2 = (dx * dx * dy * dy) / (2 * dx * dx + 2 * dy * dy);
                    grid2[i][j] = term1 * term2;
                    //pprintf("%f\t", grid2[i][j]);
                    double diff = fabs(grid2[i][j] - grid1[i][j]);
                    maxDiff = diff > maxDiff ? diff : maxDiff;
                    y += dy;
                }
                //pprintf("\n");
                x += dx;
            }
            tempGrid = grid2;
            grid2 = grid1;
            grid1 = tempGrid;
            tempGrid = NULL;
            //pprintf("MAX_DIFF: %f\n", maxDiff);
        } while (maxDiff > params->m_Tolerance);

        /* Write results */
        if (WriteGrid("laplace.dat", params, grid2))
        {
            exit(1);
        }

        /* Clean up */
        FreeGrid(grid2);
        FreeGrid(grid1);
    }
    return 0;
}

