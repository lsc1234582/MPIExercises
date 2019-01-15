#include "mpi.h"
#include "mpe.h"

#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"

#define MAX_VALUE 10

typedef enum LogEvents
{
    Bcast_Beg,
    Bcast_End,
    Compute_Beg,
    Compute_End,
    Send_Beg,
    Send_End,
    Recv_Beg,
    Recv_End,
} LogEvents;

/* Process mapping function */
int proc_map(int i, int size, int arow)
{
    size -= 1;
    int r = (int) ceil( (double) arow / (double) size);
    int proc = i / r;
    return proc + 1;
}

void print_help(void)
{
    printf("Format: main <Number of Row> <Number of Col>\n");
}

int main(int argc, char** argv)
{
    if(argc < 3)
    {
        print_help();
        exit(1);
    }

    const int AROW = atoi(argv[1]);
    const int ACOL = atoi(argv[2]);

    int size, rank;
    MPI_Status stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Initialise MPE logging
    MPE_Init_log();
    if (rank == 0)
    {
        // Describe MPE logging states
        MPE_Describe_state(Bcast_Beg, Bcast_End, "Bcast", "red:vlines3");
        MPE_Describe_state(Compute_Beg, Compute_End, "Compute", "blue:gray3");
        MPE_Describe_state(Send_Beg, Send_End, "Send", "green:light_gray");
        MPE_Describe_state(Recv_Beg, Recv_End, "Recv", "yellow:gray");

        int** a = (int**) malloc(AROW * ACOL * sizeof(int *));
        for (int i = 0; i < AROW; ++i)
        {
            a[i] = (int*) calloc(ACOL, sizeof(int));
        }
        int* b = (int*) calloc(ACOL, sizeof(int));
        int* c = (int*) calloc(AROW, sizeof(int));

        /* Generating Random Values for A & B Array*/
        srand(time(NULL));
        for (int i=0;i<AROW;i++)
        {
            for (int j=0;j<ACOL;j++)
            {
                if (i==0) b[j] = rand() % MAX_VALUE;
                a[i][j] = rand() % MAX_VALUE;
            }
        }

        /* Printing the Matrix*/
#ifdef DEBUG
        printf("Matrix A :\n");
        for (int i=0;i<AROW;i++)
        {
            for (int j=0;j<ACOL;j++)
            {
                printf("%3d ", a[i][j]);
            }
            printf("\n");
        }
        printf("nMatrix B :\n");
        for (int i=0;i<ACOL;i++)
        {
            printf("%3d ", b[i]);
        }
        printf("\n\n");
#endif

        if (size == 1)
        {
            /* Serial case */
            MPE_Log_event(Compute_Beg, 0, "compute");
            for (int i=0;i<AROW;i++)
            {
                for (int j=0; j < ACOL; ++j)
                {
                    c[i] += a[i][j] * b[j];
                }
#ifdef DEBUG
                printf("P%d : c[%d]\t= %d\n", rank, i, c[i]);
#endif
            }
            MPE_Log_event(Compute_End, 0, "computed");
        }
        else
        {
            /* Parallel case */
            /* (1) Sending B Values to other processes */
            for (int j=1;j<size;j++)
            {
                MPE_Log_event(Send_Beg, j, "send");
                MPI_Send(b, ACOL, MPI_INT, j, 99, MPI_COMM_WORLD);
                MPE_Log_event(Send_End, j, "sent");
            }

            /* (2) Sending Required A Values to specific process */
            for (int i=0;i<AROW;i++)
            {
                int processor = proc_map(i, size, AROW);
                MPE_Log_event(Send_Beg, i, "send");
                MPI_Send(a[i], ACOL, MPI_INT, processor, (100*(i+1)), MPI_COMM_WORLD);
                MPE_Log_event(Send_End, i, "sent");
            }

            /* (3) Gathering the result from other processes*/
            for (int i=0;i<AROW;i++)
            {
                int source_process = proc_map(i, size, AROW);
                MPE_Log_event(Recv_Beg, i, "recv");
                MPI_Recv(&c[i], 1, MPI_INT, source_process, i, MPI_COMM_WORLD, &stat);
                MPE_Log_event(Recv_End, i, "recvd");
#ifdef DEBUG
                printf("P%d : c[%d]\t= %d\n", rank, i, c[i]);
#endif
            }
        }

        for (int i = 0; i < AROW; ++i)
        {
            free(a[i]);
        }
        free(a);
        free(b);
        free(c);
    }
    else
    {
        int b[ACOL];

        /* (1) Each process get B Values from Master */
        MPE_Log_event(Recv_Beg, 0, "recv");
        MPI_Recv(b, ACOL, MPI_INT, 0, 99, MPI_COMM_WORLD, &stat);
        MPE_Log_event(Recv_End, 0, "recvd");

        /* (2) Get Required A Values from Master then Compute the result */
        for (int i=0;i<AROW;i++)
        {
            int processor = proc_map(i, size, AROW);
            if (rank == processor)
            {
                int buffer[ACOL];
                MPE_Log_event(Recv_Beg, i, "recv");
                MPI_Recv(buffer, ACOL, MPI_INT, 0, (100*(i+1)), MPI_COMM_WORLD, &stat);
                MPE_Log_event(Recv_End, i, "recvd");
                int sum = 0;
                MPE_Log_event(Compute_Beg, i, "compute");
                for (int j=0;j<ACOL;j++)
                {
                    sum = sum + (buffer[j] * b[j] );
                }
                MPE_Log_event(Compute_End, i, "computed");
                MPE_Log_event(Send_Beg, i, "send");
                MPI_Send(&sum, 1, MPI_INT, 0, i, MPI_COMM_WORLD);
                MPE_Log_event(Send_End, i, "sent");
            }
        }
    }
    MPE_Finish_log("pmatvec.log");
    MPI_Finalize();
    return 0;
}
