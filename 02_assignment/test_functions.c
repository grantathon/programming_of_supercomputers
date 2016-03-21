/**
 * Functions to test the data distribution and communication lists creation algorithms
 *
 * @date 22-Oct-2012, 03-Nov-2014
 * @author V. Petkov, A. Berariu
 */

#include <stdio.h>

int test_distribution(char *file_in, char *file_vtk_out, int *local_global_index,
                      int local_num_elems, double *scalars) {
    
    // global sized variables, for reading the input file
    int nintci_m, nintcf_m;  
    int nextci_m, nextcf_m;
    int **lcc_m; 
    double *bs_m, *be_m, *bn_m, *bw_m, *bh_m, *bl_m;
    double *bp_m;  
    double *su_m;  
    int points_count_m;  
    int** points_m;  
    int* elems_m;  
    int i;

    // read the entire file
    int f_status = read_binary_geo( file_in, &nintci_m, &nintcf_m, &nextci_m,
            &nextcf_m, &lcc_m, &bs_m, &be_m, &bn_m, &bw_m, &bl_m, &bh_m, &bp_m,
            &su_m, &points_count_m, &points_m, &elems_m );
    if ( f_status != 0 ) {
        printf( "Error in reading input file \n" );
        return -1;
    }

    // allocate distribution vector
    double *distr;
    if ( ( distr = (double *) malloc( ( nintcf_m + 1 ) * sizeof(double) ) )
            == NULL ) {
        printf( "malloc failed to allocate distr array" );
        return -1;
    }
    for ( i = nintci_m; i < ( nintcf_m + 1 ); i++ ) {
        distr[i] = 0.0;
    }

    // copy the local values using the generated map
    for ( i = 0; i < local_num_elems; i++ ) {
        distr[local_global_index[i]] = scalars[i];
    }
    // write vtk file
    vtk_write_unstr_grid_header( file_in, file_vtk_out, nintci_m, nintcf_m,
            points_count_m, points_m, elems_m );
    vtk_append_double( file_vtk_out, "SCALARS", nintci_m, nintcf_m, distr );
    printf( "Distribution VTK file succesfully generated! \n" );

    // free the allocated memory
    free( su_m );
    free( bp_m );
    free( bh_m );
    free( bl_m );
    free( bw_m );
    free( bn_m );
    free( be_m );
    free( bs_m );
    free( elems_m );

    for ( i = 0; i < nintcf_m + 1; i++ ) {
        free( lcc_m[i] );
    }
    free( lcc_m );

    for ( i = 0; i < points_count_m; i++ ) {
        free( points_m[i] );
    }
    free( points_m );
    free( distr );

    return 0;
}

/* Write statistics to pstats.dat 
 * @param input_key: 1 - tjunc
 					 2 - drall
 					 3 - pent
 					 4 - cojack
 * @param part_key: 1 - classic
 					2 - dual
 					3 - nodal
 * @param read_key: 1 - oneread
 				    2 - allread
 * @param my_rank: current process rank
 * @param time_usec: execution time in microseconds
 * @return
 */

int write_pstats_exectime( int input_key,
				           int part_key,
				           int read_key,
				           int my_rank,
				           long long time_usec ){

	// append to existing file
    char pstatsname[20];
    sprintf(pstatsname, "pstats_%d.dat", my_rank);
    FILE *fp = fopen( pstatsname, "a" );
    if( fp == NULL ){
        printf( "Error opening file %s for writing\n", pstatsname );
        return -1;
    }


    fprintf( fp, "EXECTIME %d %d %d %d %lld\n", input_key, part_key, read_key, my_rank, time_usec );
    fclose( fp );

    return 0;
}

/* Write statistics to pstats.dat 
 * @param input_key: 1 - tjunc
 					 2 - drall
 					 3 - pent
 					 4 - cojack
 * @param part_key: 1 - classic
 					2 - dual
 					3 - nodal
 * @param my_rank: current process rank
 * @param local_intc: number of local internal cells
 * @param local_extc: number of local external cells
 * @return
 */

int write_pstats_partition( int input_key,
				            int part_key,
				            int my_rank,
				            int local_intc,
				            int local_extc ){

	// append to existing file
    char pstatsname[20];
    sprintf(pstatsname, "pstats_%d.dat", my_rank);
    FILE *fp = fopen( pstatsname, "a" );
    if( fp == NULL ){
        printf( "Error opening file %s for writing\n", pstatsname );
        return -1;
    }


    fprintf( fp, "PARTITION %d %d %d %d %d\n", input_key, part_key, my_rank, local_intc, local_extc );
    fclose( fp );

    return 0;
}

/* Write statistics to pstats.dat 
 * @param input_key: 1 - tjunc
                     2 - drall
                     3 - pent
                     4 - cojack
 * @param part_key: 1 - classic
                    2 - dual
                    3 - nodal
 * @param my_rank: current process rank
 * @param nprocs: the total number of processes
 * @param nghb_cnt: the total number of neighbours
 * @param nghb_idx: the index of the neighbour to output data for
 * @param send_cnt: the number of cells to be sent to each neighbour
 * @param send_lst: the indexes of the cells to be sent to each neighbour
 * @param recv_cnt: the number of cells to be received from each neighbour
 * @param recv_lst: the indexes of the cells to be received from each neighbour
 * @return
 */

int write_pstats_communication( int input_key,
                            int part_key,
                            int my_rank,
                            int nprocs,
                            int nghb_cnt,
                            int nghb_idx,
                            int* send_cnt,
                            int** send_lst,
                            int* recv_cnt,
                            int** recv_lst ){

    // append to existing file
    char pstatsname[20];
    sprintf(pstatsname, "pstats_%d.dat", my_rank);
    FILE *fp = fopen( pstatsname, "a" );
    if( fp == NULL ){
        printf( "Error opening file %s for writing\n", pstatsname );
        return -1;
    }


    fprintf( fp, "COMMUNICATION %d %d %d %d %d %d %d %d %d %d %d %d\n", 
        input_key, part_key, my_rank, nprocs, nghb_cnt, nghb_idx, 
        send_lst[nghb_idx][0], send_lst[nghb_idx][(send_cnt[nghb_idx]-1)/2], send_lst[nghb_idx][send_cnt[nghb_idx]-1],
        recv_lst[nghb_idx][0], recv_lst[nghb_idx][(recv_cnt[nghb_idx]-1)/2], recv_lst[nghb_idx][recv_cnt[nghb_idx]-1]);
    fclose( fp );

    return 0;
}

