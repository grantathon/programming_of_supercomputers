/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include "mpi.h"
#include "util_write_files.h"

void finalization(char* file_in, int nprocs, int myrank, int total_iters, double residual_ratio,
                  int nintci, int nintcf, double* var, int* local_global_index, int* global_local_index)
{
    // int global_nintci, global_nintcf;
    int i, j, status, global_int_count, local_int_count;
    double *unordered_global_vars, *ordered_global_vars;
    int *int_counts, *int_disps, *local_global_index_all;
    char file_out[100];

    sprintf(file_out, "%s.summary.out", file_in);

    // Perfrom inter-process communication if necessary
    if(nprocs > 1)
    {
        local_int_count = nintcf - nintci + 1;
        int_counts = (int *)malloc(nprocs * sizeof(int));
        int_disps = (int *)malloc(nprocs * sizeof(int));

        // Communicate neighboring procs' internal cell counts
        MPI_Gather(&local_int_count, 1, MPI_INT, int_counts, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if(myrank == 0)
        {
            global_int_count = 0;
            for(i = 0; i < nprocs; i++)
            {
                int_disps[i] = global_int_count;
                global_int_count += int_counts[i];
            }

            local_global_index_all = (int *)malloc(global_int_count * sizeof(int));
            unordered_global_vars = (double *)malloc(global_int_count * sizeof(double));
        }

        // Communicate neighboring procs' local-to-global index mappings
        MPI_Gatherv(local_global_index, (nintcf + 1), MPI_INT, local_global_index_all, int_counts, int_disps, MPI_INT, 0, MPI_COMM_WORLD);

        // Communicate neighboring procs' var elements
        MPI_Gatherv(var, (nintcf + 1), MPI_DOUBLE, unordered_global_vars, int_counts, int_disps, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        if(myrank == 0)
        {
            ordered_global_vars = (double *)malloc(global_int_count * sizeof(double));

            // Input vars from master rank
            for(j = 0; j < int_counts[0]; j++)
            {
                ordered_global_vars[local_global_index_all[j]] = unordered_global_vars[j];
            }

            // Input vars from neighbor ranks
            for(i = 1; i < nprocs; i++)
            {
                for(j = 0; j < int_counts[i]; j++)
                {
                    //ordered_global_vars[local_global_index_all[j + int_counts[i - 1]]] = unordered_global_vars[j + int_counts[i - 1]];
                    ordered_global_vars[local_global_index_all[j + int_disps[i]]] = unordered_global_vars[j + int_disps[i]];
                }
            }

            status = store_simulation_stats(file_in, file_out, 0, (global_int_count - 1), ordered_global_vars, total_iters, residual_ratio);

            if ( status != 0 ) fprintf(stderr, "Error when trying to write to file %s\n", file_out);

            // Deallocate memory
            free(ordered_global_vars);
            free(unordered_global_vars);
            free(local_global_index_all);
        }

        // Deallocate memory
        free(int_counts);
        free(int_disps);
    }
    else
    {
        status = store_simulation_stats(file_in, file_out, nintci, nintcf, var, total_iters, residual_ratio);

        if ( status != 0 ) fprintf(stderr, "Error when trying to write to file %s\n", file_out);
    }
}

