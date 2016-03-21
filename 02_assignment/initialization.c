/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012, 13-Nov-2014
 * @author V. Petkov, A. Berariu
 */

#include <stdlib.h>
#include "mpi.h"

#include "util_read_files.h"
#include "initialization.h"
#include "test_functions.h"
 #include "data_distribution.h"

int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank,
                   int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index, int** global_local_index,
                   int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, 
                   int **recv_cnt, int*** recv_lst)
{
    int i = 0;

    // locall cell counts
    int f_status, local_int_count, local_ext_count;

    // variables for writing stats
    int input_key, part_key, read_key;

    // read-in the input file
    f_status = read_data(   file_in, part_type, read_type, nintci, nintcf, nextci, nextcf,
                            lcc, bs, be, bn, bw, bl, bh, bp, su, points_count, points, elems,
                            local_global_index, global_local_index, nprocs, myrank, &local_int_count,
                            &local_ext_count, nghb_cnt, nghb_to_rank, send_cnt, send_lst, 
                            recv_cnt, recv_lst);

    if ( f_status != 0 ) return f_status;

    *var = (double*) calloc(sizeof(double), (*nextcf + 1));
    *cgup = (double*) calloc(sizeof(double), (*nextcf + 1));
    *cnorm = (double*) calloc(sizeof(double), (*nintcf + 1));

    // initialize the arrays
    for ( i = 0; i <= 10; i++ ) {
        (*cnorm)[i] = 1.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ ) {
        (*var)[i] = 0.0;
    }

    for ( i = (*nintci); i <= (*nintcf); i++ ) {
        (*cgup)[i] = 1.0 / ((*bp)[i]);
    }

    for ( i = (*nextci); i <= (*nextcf); i++ ) {
        (*var)[i] = 0.0;
        (*cgup)[i] = 0.0;
        (*bs)[i] = 0.0;
        (*be)[i] = 0.0;
        (*bn)[i] = 0.0;
        (*bw)[i] = 0.0;
        (*bh)[i] = 0.0;
        (*bl)[i] = 0.0;
    }

    return 0;
}

