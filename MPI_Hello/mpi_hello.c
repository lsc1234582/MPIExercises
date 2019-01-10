/******************************************************************************
 * * FILE: mpi_hello.c
 * * DESCRIPTION:
 * *   MPI tutorial example code: Simple hello world program
 * * AUTHOR: Blaise Barney
 * * LAST REVISED: 03/05/10
 * ******************************************************************************/
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#define  MASTER     0

int main (int argc, char *argv[])
{
    int   numtasks, taskid, len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    // Initialising the MPI execution environment
    MPI_Init(&argc, &argv);
    // Get the size of a group (I.e. the number of processes within a group), identified by the communicator
    // "MPI_COMM_WORLD".
    // A communicator is an identifier that indentifies a group of processes and a communication context.
    // A communication context, or simply a context, is similar to tags in that it is used to match messages, but
    // differs in that it's allocated by the system instead of the user.
    // A tag is a user-allocated ID used by both the sender and the receiver to match messages, similar to
    // the notion of "channels".
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    // Get the rank of a process within the group associated with the communicator "MPI_COMM_WORLD".
    // The rank of a process is its ID within its group, and it ranges from 0 to group size - 1.
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    // Get the name of the processor.
    MPI_Get_processor_name(hostname, &len);
    printf ("Hello from task %d on %s!\n", taskid, hostname);
    if (taskid == MASTER)
    {
        // Print a message in process with rank MASTER.
        printf("MASTER: Number of MPI tasks is: %d\n",numtasks);
    }
    // Terminates MPI execution environment.
    MPI_Finalize();

}
