/**
 * Main GCCG program
 *
 * @author E. Xue, V. Petkov, A. Berariu
 * @date 22-May-2009, 22-Oct-2012, 13-Nov-2014
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "initialization.h"
#include "compute_solution.h"
#include "finalization.h"
#include "test_functions.h"

int main(int argc, char *argv[])
{
    int my_rank, num_procs, i;
    const int max_iters = 10000;    /// maximum number of iteration to perform

    /** Simulation parameters parsed from the input datasets */
    int nintci, nintcf;    /// internal cells start and end index

    /// external cells start and end index. The external cells are only ghost cells.
    /// They are accessed only through internal cells
    int nextci, nextcf;
    int **lcc;    /// link cell-to-cell array - stores neighboring information

    /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
    double *bs, *be, *bn, *bw, *bh, *bl;
    double *bp;    /// Pole coefficient
    double *su;    /// Source values
    double residual_ratio;    /// the ratio between the reference and the current residual
    double *var;    /// the variation vector -> keeps the result in the end

    /** Additional vectors required for the computation */
    double *cgup, *oc, *cnorm;

    /** Geometry data */
    int points_count;    /// total number of points that define the geometry
    int** points;    /// coordinates of the points that define the cells - size [points_cnt][3]
    int* elems;    /// definition of the cells using their nodes (points) - each cell has 8 points

    /** Mapping between local and remote cell indices */
    int* local_global_index;    /// local to global index mapping
    int* global_local_index;    /// global to local index mapping
  
    /** Lists for neighbouring information */
    int nghb_cnt = 0;    /// total number of neighbors of the current process
    int *nghb_to_rank;  /// mapping of the neighbour index to the corresponding process rank
    int *send_cnt;    /// number of cells to be sent to each neighbour (size: nghb_cnt)
    int **send_lst;    /// lists of cells to be sent to each neighbour (size: nghb_cnt x send_cnt[*])
    int *recv_cnt;    /// number of cells to be received from each neighbour (size: nghb_cnt)
    int **recv_lst;    /// lists of cells to be received from each neighbour (size: nghb_cnt x recv_cnt[*])

    // MPI initialization functions
    MPI_Init(&argc, &argv);    /// Start MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /// get current process id
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);    /// get number of processes

    /** process call arguments **/
    if ( argc < 4 ) {
        fprintf(stderr, "Usage: ./gccg <input_file> <partition_type> <algorithm_type>\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    char *file_in = argv[1];
    char *part_type = argv[2];
    if ( strcmp( part_type, "classic" ) && strcmp( part_type, "dual" )
            && strcmp( part_type, "nodal" ) ) {
        printf(
                " Wrong partition type selected. Valid values are classic, nodal and dual \n" );
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    char *read_type = argv[3];
    if ( strcmp( read_type, "oneread" ) && strcmp( read_type, "allread" ) ) {
        printf(
                " Wrong read-in algorithm selected. Valid values are oneread and allread. \n" );
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    /********** START INITIALIZATION **********/
    
    // read-in the input file
    int init_status = initialization(file_in, part_type, read_type, num_procs, my_rank,
                                     &nintci, &nintcf, &nextci, &nextcf, 
                                     &lcc, &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, 
                                     &points_count, &points, &elems, &var, &cgup, &oc, &cnorm, 
                                     &local_global_index, &global_local_index,
                                     &nghb_cnt, &nghb_to_rank, 
                                     &send_cnt, &send_lst, &recv_cnt, &recv_lst);

    if ( init_status != 0 ) {
        fprintf(stderr, "Failed to initialize data!\n");
        MPI_Abort(MPI_COMM_WORLD, my_rank);
    }

    /********** END INITIALIZATION **********/

    /********** START COMPUTATIONAL LOOP **********/

    // int total_iters = compute_solution_serial(num_procs, my_rank, max_iters, nintci, nintcf, nextcf, 
    //                                           lcc, bp, bs, bw, bl, bn, be, bh,
    //                                           cnorm, var, su, cgup, &residual_ratio,
    //                                           local_global_index, global_local_index, nghb_cnt, 
    //                                           nghb_to_rank, send_cnt, send_lst, recv_cnt, recv_lst);

    int total_iters = compute_solution_parallel(num_procs, my_rank, max_iters, nintci, nintcf, nextcf, 
                                                lcc, bp, bs, bw, bl, bn, be, bh,
                                                cnorm, var, su, cgup, &residual_ratio,
                                                local_global_index, global_local_index, nghb_cnt, 
                                                nghb_to_rank, send_cnt, send_lst, recv_cnt, recv_lst);

    /********** END COMPUTATIONAL LOOP **********/

    /********** START FINALIZATION **********/

    finalization(file_in, num_procs, my_rank, total_iters, residual_ratio, nintci, nintcf, var, local_global_index, global_local_index);

    /********** END FINALIZATION **********/

    // cleanup allocated memory
    free(cnorm);
    free(var);
    free(cgup);
    free(su);
    free(bp);
    free(bh);
    free(bl);
    free(bw);
    free(bn);
    free(be);
    free(bs);
    free(elems);

    for ( i = 0; i < nintcf + 1; i++ ) {
        free(lcc[i]);
    }
    free(lcc);

    for ( i = 0; i < points_count; i++ ) {
        free(points[i]);
    }
    free(points);

    free(nghb_to_rank);
    
    free(send_cnt);
    for ( i = 0; i < nghb_cnt; i++ ){
        free(send_lst[i]);
    }
    free(send_lst);

    free(recv_cnt);
    for ( i = 0; i < nghb_cnt; i++ ){
        free(recv_lst[i]);
    }
    free(recv_lst);

    free(local_global_index);
    free(global_local_index);

    MPI_Finalize();    /// cleanup MPI

    return 0;
}

