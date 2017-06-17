#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define MASTER 0
#define MAX_STEPS 1000
#define TAM 450

int mod (int a, int b)
{
	int ret = a % b;

	if (ret < 0)
		ret += b;
	return ret;
}

void processState (int *old, int *new, int lin, int col)
{
	int alive;
	int l1 = (TAM * mod (lin - 1, TAM)), l2 = (TAM * mod (lin, TAM)), l3 = (TAM * mod (lin + 1, TAM));
	int c1 = mod (col - 1, TAM), c2 =  mod (col, TAM), c3 = mod (col + 1, TAM);

	alive = old[l1 + c1] + old[l1 + c2] + old[l1 + c3] + old[l2 + c1] + old[l2 + c3] + old[l3 + c1] + old[l3 + c2] + old[l3 + c3];

	if (alive < 2)
		new[TAM * lin + col] = 0;
	else
		if (alive == 2 || alive == 3)
			new[TAM * lin + col] = 1;
		else
			new[TAM * lin + col] = 0;
}

void inicializaMatriz (int *old)
{        
	int i, j;

	srand (0);
	for (i = 0; i < TAM; ++i)
		for (j = 0; j < TAM; ++j)
			old[i * TAM + j] = rand () % 2;
}

void printMatriz (int *old)
{
        int i, j;

        for (i = 0; i < TAM; ++i)
        {
                for (j = 0; j < TAM; ++j)
                {
                        if (old[i * TAM + j] == 1)
                                printf ("%c ", 43);
                        else
                                printf ("0 ");
                }
                printf ("\n");
        }
        printf ("\n\n");
}

int main (int argc, char **argv)
{
	double thread_chunk;
	int i, j, STEPS;
	int offset, chunk;
	int *aux, *new, *old;
	int my_rank, n_procs;
	MPI_Status status;

	MPI_Init (&argc, &argv);
	MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size (MPI_COMM_WORLD, &n_procs);
	thread_chunk = TAM * 1.0 / n_procs * 1.0;

	if (my_rank == MASTER)
	{
		new = malloc (TAM * TAM * sizeof (int));
		old = malloc (TAM * TAM * sizeof (int));
		inicializaMatriz (old);

		for (i = 1; i < n_procs; ++i)
		{
			chunk = (int) floor ((i + 1) * thread_chunk) - (int) floor (i * thread_chunk);
			offset = (int) floor (i * thread_chunk);
			MPI_Send (&old[TAM * offset], TAM * chunk, MPI_INT, i, 1, MPI_COMM_WORLD);
		}

		for (STEPS = 0; STEPS < MAX_STEPS; ++STEPS)
		{
			for (i = 1; i < n_procs; ++i)
			{
				chunk = (int) floor ((i + 1) * thread_chunk) - (int) floor (i * thread_chunk);
				offset = (int) floor (i * thread_chunk);

				if (i != n_procs - 1)
				{
					MPI_Send (&old[TAM * offset - TAM], TAM, MPI_INT, i, 1, MPI_COMM_WORLD);
					MPI_Send (&old[TAM * offset + TAM * chunk], TAM, MPI_INT, i, 2, MPI_COMM_WORLD);
				}
				else
				{
					MPI_Send (&old[TAM * offset - TAM], TAM, MPI_INT, i, 1, MPI_COMM_WORLD);
					MPI_Send (&old[0], TAM, MPI_INT, i, 2, MPI_COMM_WORLD);
				}
			}

			for (i = 0; i < (int) floor (thread_chunk); ++i)
				for (j = 0; j < TAM; ++j)
					processState (old, new, i, j);

			for (i = 1; i < n_procs; ++i)
			{
				chunk = (int) floor ((i + 1) * thread_chunk) - (int) floor (i * thread_chunk);
				offset = (int) floor (i * thread_chunk);
				MPI_Recv (&new[TAM * offset], TAM, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
				MPI_Recv (&new[TAM * offset + TAM * chunk - TAM], TAM, MPI_INT, i, 2, MPI_COMM_WORLD, &status);
			}

			aux = old;
			old = new;
			new = aux;
		}

		for (i = 1; i < n_procs; ++i)
		{
			chunk = (int) floor ((i + 1) * thread_chunk) - (int) floor (i * thread_chunk);
			offset = (int) floor (i * thread_chunk);
			MPI_Recv (&old[TAM * offset], TAM * chunk, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
		}

		printMatriz (old);
	}
	else
	{
		STEPS = 0;
		chunk = (int) floor ((my_rank + 1) * thread_chunk) - (int) floor (my_rank * thread_chunk);
		new = malloc (TAM * (chunk + 2) * sizeof (int));
		old = malloc (TAM * (chunk + 2) * sizeof (int));
		MPI_Recv (&old[TAM], TAM * chunk, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
		do
		{
			++STEPS;
			MPI_Recv (&old[0], TAM, MPI_INT, MASTER, 1, MPI_COMM_WORLD, &status);
			MPI_Recv (&old[TAM * chunk + TAM], TAM, MPI_INT, MASTER, 2, MPI_COMM_WORLD, &status);

			for (i = 1; i < chunk + 1; ++i)
				for (j = 0; j < TAM; ++j)
					processState (old, new, i, j);

			MPI_Send (&new[TAM], TAM, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
			MPI_Send (&new[TAM * chunk], TAM, MPI_INT, MASTER, 2, MPI_COMM_WORLD);

			aux = old;
			old = new;
			new = aux;
		} while (STEPS < MAX_STEPS);
		MPI_Send (&old[TAM], TAM * chunk, MPI_INT, MASTER, 1, MPI_COMM_WORLD);
	}
	MPI_Finalize ();
}
