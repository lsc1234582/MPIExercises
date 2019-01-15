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
            int sentItemsPlusOne = 1;
            /* Parallel case */
            /* (1) Broadcasting B Values to other processes */
            MPE_Log_event(Bcast_Beg, 0, "broadcast");
            MPI_Bcast(b, ACOL, MPI_INT, 0, MPI_COMM_WORLD);
            MPE_Log_event(Bcast_End, 0, "broadcasted");

            /* (2) Sending Required A Values to specific process */
            for (int i = 1; i < size; ++i)
            {
                MPE_Log_event(Send_Beg, i, "send");
                MPI_Send(a[i], ACOL, MPI_INT, i, sentItemsPlusOne, MPI_COMM_WORLD);
                MPE_Log_event(Send_End, i, "sent");
                ++sentItemsPlusOne;
                if(sentItemsPlusOne == AROW + 1)
                {
                    break;
                }
            }

            /* (3) Gathering the result from other processes*/
            for (int i=0; i<AROW; i++)
            {
                int resElem;
                MPE_Log_event(Recv_Beg, i, "recv");
                MPI_Recv(&resElem, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
                MPE_Log_event(Recv_End, i, "recvd");
                int resIndexPlusOne = stat.MPI_TAG;
                c[resIndexPlusOne - 1] = resElem;
#ifdef DEBUG
                printf("P%d : c[%d]\t= %d\n", rank, resIndexPlusOne - 1, c[resIndexPlusOne - 1]);
#endif
                /* If there're more items, dispatch them */
                if (sentItemsPlusOne < AROW + 1)
                {
                    for (int i = 1; i < size; ++i)
                    {
                        MPE_Log_event(Send_Beg, i, "send");
                        MPI_Send(a[i], ACOL, MPI_INT, i, sentItemsPlusOne, MPI_COMM_WORLD);
                        MPE_Log_event(Send_End, i, "sent");
                        ++sentItemsPlusOne;
                        if(sentItemsPlusOne == AROW + 1)
                        {
                            break;
                        }
                    }
                }
            }

            /* (4) Inform all processes jobs' done */
            for (int i = 1; i < size; ++i)
            {
                // message jobsDone is not used
                int jobsDone = 0;
                MPE_Log_event(Send_Beg, i, "send");
                MPI_Send(&jobsDone, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPE_Log_event(Send_End, i, "sent");
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
        int receivedItems = 0;

        /* (1) Receiving B Values from Master broadcast */
        MPE_Log_event(Bcast_Beg, 0, "broadcast");
        MPI_Bcast(b, ACOL, MPI_INT, 0, MPI_COMM_WORLD);
        MPE_Log_event(Bcast_End, 0, "broadcasted");

        /* (2) Get Required A Values from Master then Compute the result */
        while (1)
        {
            int buffer[ACOL];
            MPE_Log_event(Recv_Beg, receivedItems, "recv");
            MPI_Recv(buffer, ACOL, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
            MPE_Log_event(Recv_End, receivedItems, "recvd");
            int resIndexPlusOne = stat.MPI_TAG;
            /* If received tag is 0, then stop */
            if (resIndexPlusOne == 0)
            {
                break;
            }
            int sum = 0;
            MPE_Log_event(Compute_Beg, receivedItems, "compute");
            for (int j=0;j<ACOL;j++)
            {
                sum = sum + (buffer[j] * b[j] );
            }
            MPE_Log_event(Compute_End, receivedItems, "computed");

            MPE_Log_event(Send_Beg, receivedItems, "send");
            MPI_Send(&sum, 1, MPI_INT, 0, resIndexPlusOne, MPI_COMM_WORLD);
            MPE_Log_event(Send_End, receivedItems, "sent");

            ++receivedItems;
        }
    }
    MPE_Finish_log("pmatvec.log");
    MPI_Finalize();
    return 0;
}
