/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2004 by University of Chicago.
 *      See COPYRIGHT in top-level directory.
 */

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "mlife2d.h"

int MLIFE_ExchangeInitPt2ptVector( MLIFEPatchDesc *patch,
                 int **m1, int **m2, void *privateData )
{
    *(void **)privateData = 0;
    return 0;
}

int MLIFE_ExchangeEndPt2ptVector( void *privateData )
{
    return 0;
}

int MLIFE_ExchangePt2ptVector( MLIFEPatchDesc *patch, int **matrix,
               MLIFETiming *timedata, void *privateData )
{
    MPI_Request reqs[4];
    MPI_Comm    comm = patch->comm;
    int LRows = patch->lni;
    int LCols = patch->lnj;

    /* Send and receive boundary information */
    /*
     *        |<- lnj            ->|      
     * ------------------------------------ -
     * || lr1 | ls1 | eo ... | rs1 | rr1 || 
     * || lr2 | ls2 |  ...   | rs2 | rr2 || lni
     * || ..  | ..  |  ...   | ..  | ..  ||
     * || lrm | lsm |  ...   | rsm | rrm ||
     * ------------------------------------ -
     */
    MPI_Datatype VertEdgeType;
    MPI_Type_vector(LRows, 1, LCols + 2, MPI_INT, &VertEdgeType);
    MPI_Type_commit(&VertEdgeType);

    /* first, move the left, right edges */
    MPI_Isend(&(matrix[1][1]), 1, VertEdgeType, patch->left, 0, comm, reqs);
    MPI_Irecv(&(matrix[1][0]), 1, VertEdgeType, patch->left, 0, comm, reqs+1);
    MPI_Isend(&(matrix[1][LCols]), 1, VertEdgeType, patch->right, 0, comm, reqs+2);
    MPI_Irecv(&(matrix[1][LCols + 1]), 1, VertEdgeType, patch->right, 0, comm, reqs+3);
    /* We need to wait on these for the trick that we use to move
       the diagonal terms to work */
    MPI_Waitall( 4, reqs, MPI_STATUSES_IGNORE );
    /* move the top, bottom edges (including diagonals) */
    MPI_Isend(&matrix[1][0], LCols+2, MPI_INT,
          patch->up, 0, comm, reqs);
    MPI_Irecv(&matrix[0][0], LCols+2, MPI_INT,
          patch->up, 0, comm, reqs+1);
    MPI_Isend(&matrix[LRows][0], LCols+2, MPI_INT,
          patch->down, 0, comm, reqs+2);
    MPI_Irecv(&matrix[LRows+1][0], LCols+2, MPI_INT,
          patch->down, 0, comm, reqs+3);

    MPI_Waitall(4, reqs, MPI_STATUSES_IGNORE);

    return MPI_SUCCESS;
}
