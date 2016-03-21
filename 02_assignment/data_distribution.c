#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

#include "data_distribution.h"
#include "metis.h"

#define TRUE 1
#define FALSE 0

int read_data(  char *file_name, char *part_type, char *read_type, int *nintci, int *nintcf, int *nextci,
                int *nextcf, int ***lcc, double **bs, double **be, double **bn, double **bw, double **bl,
                double **bh,double **bp, double **su, int *points_count, int ***points, int **elems,
                int **local_global_index, int **global_local_index, int nprocs, int myrank, int *local_int_count_out,
                int *local_ext_count_out, int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, 
                int **recv_cnt, int*** recv_lst)
{
    int f_status = 0;

    // Determine the read type
    if(!strcmp(read_type, "oneread"))
    {
        f_status = one_proc_read(   file_name, part_type, nintci, nintcf, nextci, nextcf, lcc, bs, be,
                                    bn, bw, bl, bh, bp, su, points_count, points, elems, local_global_index,
                                    global_local_index, nprocs, myrank, local_int_count_out, local_ext_count_out,
                                    nghb_cnt, nghb_to_rank, send_cnt, send_lst, recv_cnt, recv_lst);
    }
    else if(!strcmp(read_type, "allread"))
    {
        f_status = all_proc_read(   file_name, part_type, nintci, nintcf, nextci, nextcf, lcc, bs, be,
                                    bn, bw, bl, bh, bp, su, points_count, points, elems, local_global_index,
                                    global_local_index, nprocs, myrank, local_int_count_out, local_ext_count_out,
                                    nghb_cnt, nghb_to_rank, send_cnt, send_lst, recv_cnt, recv_lst);
    }
    else
    {
        if(myrank == 0)
        {
            printf("ERROR: Unsupported input file reading type: %s\n", read_type);  
        }

        return 1;
    }


    if ( f_status != 0 ) return f_status;

    return 0;
}

int all_proc_read(  char *file_name, char *part_type, int *nintci, int *nintcf, int *nextci, int *nextcf,
                    int ***lcc, double **bs, double **be, double **bn, double **bw, double **bl,
                    double **bh,double **bp, double **su, int *points_count, int ***points, int **elems,
                    int **local_global_index, int **global_local_index, int nprocs, int myrank, int *local_int_count_out,
                    int *local_ext_count_out, int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, 
                    int **recv_cnt, int*** recv_lst)
{
    int f_status, i, j, k, cnt, is_new_rank;
    int *elem_part, *epart;
    int global_int_start_idx, global_double_start_idx, points_start_idx;
    int global_nintci, global_nintcf, global_nextci, global_nextcf;
    int **global_lcc;
    double *global_bs, *global_be, *global_bn, *global_bw, *global_bl, *global_bh, *global_bp,
        *global_su;
    int global_int_count, global_ext_count, local_int_count, std_int_count, local_ext_count, rem_ints,
        rem_exts, points_rem;

    // Read all data from input file
    f_status = read_binary_geo(file_name, &global_nintci, &global_nintcf, &global_nextci, &global_nextcf,
                            &global_lcc, &global_bs, &global_be, &global_bn, &global_bw, &global_bl, &global_bh,
                            &global_bp, &global_su, points_count, points, elems);

    if ( f_status != 0 ) return f_status;

    // Determine number of ints and doubles to distribute
    global_int_count = global_nintcf - global_nintci + 1;
    global_ext_count = global_nextcf - global_nextci + 1;

    epart = (int *)malloc(global_int_count * sizeof(int));

    if(!strcmp(part_type, "classic"))
    {
        // Index partitioning parameters
        local_int_count = global_int_count / nprocs;
        std_int_count = local_int_count;
        rem_ints = global_int_count % nprocs;

        // Determine subdomain starting indices to assign to current processor
        global_int_start_idx = myrank*local_int_count + global_nintci;

        // Add remaining ints and doubles to last processes (if necessary)
        if(myrank == nprocs - 1)
        {
            local_int_count += rem_ints;
        }

        // Determine the global index for each local cell
        elem_part = (int *)malloc(local_int_count * sizeof(int));
        for(i = 0; i < local_int_count; i++)
        {
            elem_part[i] = global_int_start_idx + i;
        }

        // Specify the element partition for each processor
        cnt = std_int_count;
        for(i = 0; i < nprocs; i++)
        {
            // Change the count for the last process to get remainder elements
            if(i == nprocs - 1)
            {
                cnt = std_int_count + rem_ints;
            }

            for(j = 0; j < cnt; j++)
            {
                epart[i*std_int_count + j] = i;
            }
        }
    }
    else
    {
        idx_t *epart_metis = (idx_t *)malloc( (idx_t)global_int_count * sizeof(idx_t));

        f_status = metis_partition(part_type, global_int_count, *points_count, nprocs, *elems, &epart_metis);

        if ( f_status != 0 ) return f_status;

        elem_part = (int *)malloc(global_int_count * sizeof(int));

        for (i = 0; i < global_int_count; i++)
        {
            epart[i] = epart_metis[i];
        }

        // Convert METIS output to local structure
        local_int_count = 0;
        for(i = 0; i < global_int_count; i++)
        {
            if(epart[i] == myrank)
            {
                elem_part[local_int_count++] = i;
            }
        }

        free(epart_metis);
    }

    // Map local to global indices
    map_local_global_indices(local_global_index, elem_part, local_int_count);

    // Allocate memory and retrieve subdomain data
    *lcc = (int**) malloc(local_int_count * sizeof(int*));
    for(i = 0; i < local_int_count; i++)
    {
        (*lcc)[i] = (int *) malloc(6 * sizeof(int));
    }

    // Setup lcc and count number of external cells
    local_ext_count = 0;
    for(i = 0; i < local_int_count; i++)
    {
        for (j = 0; j < 6; j++)
        {
            (*lcc)[i][j] = global_lcc[(*local_global_index)[i]][j];
            
            if((*lcc)[i][j] > global_nintcf) local_ext_count++;
        }
    }

    find_neighbor_procs(*lcc, epart, nghb_cnt, nghb_to_rank, global_nintcf, nprocs, local_int_count, myrank);

    // find the send count for each neighbour
    recvSendCountList(*lcc, epart, (*nghb_cnt), *nghb_to_rank, send_cnt, recv_cnt, send_lst, recv_lst, 
                        local_int_count, myrank, global_nintcf, *local_global_index, *global_local_index);

    // Set local index boundaries
    *nintci = 0;
    *nintcf = local_int_count - 1;
    *nextci = local_int_count;
    *nextcf = local_ext_count + local_int_count - 1;

    // sort
    sortListForAll((*nghb_cnt), *recv_cnt, recv_lst, myrank, *global_local_index);
    sortListForAll((*nghb_cnt), *send_cnt, send_lst, myrank, *global_local_index);

    // Map global to local indices
    map_global_local_indices(global_local_index, epart, *lcc, *recv_cnt, *recv_lst, *nghb_cnt, global_int_count, global_ext_count, local_int_count, global_nextci, myrank);

    *bs = (double *) malloc((local_ext_count + local_int_count) * sizeof(double));
    *be = (double *) malloc((local_ext_count + local_int_count) * sizeof(double));
    *bn = (double *) malloc((local_ext_count + local_int_count) * sizeof(double));
    *bw = (double *) malloc((local_ext_count + local_int_count) * sizeof(double));
    *bl = (double *) malloc((local_ext_count + local_int_count) * sizeof(double));
    *bh = (double *) malloc((local_ext_count + local_int_count) * sizeof(double));
    *bp = (double *) malloc((local_ext_count + local_int_count) * sizeof(double));
    *su = (double *) malloc((local_ext_count + local_int_count) * sizeof(double));

    for(i = 0; i < local_int_count; i++)
    {
        (*bs)[i] = global_bs[(*local_global_index)[i]];
        (*be)[i] = global_be[(*local_global_index)[i]];
        (*bn)[i] = global_bn[(*local_global_index)[i]];
        (*bw)[i] = global_bw[(*local_global_index)[i]];
        (*bl)[i] = global_bl[(*local_global_index)[i]];
        (*bh)[i] = global_bh[(*local_global_index)[i]];
        (*bp)[i] = global_bp[(*local_global_index)[i]];
        (*su)[i] = global_su[(*local_global_index)[i]];
    }

    if ( f_status != 0 ) return f_status;


    (*local_int_count_out) = local_int_count;
    (*local_ext_count_out) = local_ext_count;

    // Free all temporary data structures
    free(global_bs);
    free(global_be);
    free(global_bn);
    free(global_bw);
    free(global_bl);
    free(global_bh);
    free(global_bp);
    free(global_su);
    free(elem_part);
    free(epart);

    for(i = global_nintci; i < global_nintcf + 1; i++)
    {
        free(global_lcc[i]);
    }
    free(global_lcc);

    return 0;
}

int one_proc_read(  char *file_name, char *part_type, int *nintci, int *nintcf, int *nextci, int *nextcf,
                    int ***lcc, double **bs, double **be, double **bn, double **bw, double **bl,
                    double **bh, double **bp, double **su, int *points_count, int ***points, int **elems,
                    int **local_global_index, int **global_local_index, int nprocs, int myrank, int *local_int_count_out,
                    int *local_ext_count_out, int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, 
                    int **recv_cnt, int*** recv_lst)
{
    int f_status, i, j, cnt;
    int local_int_count_min, local_ext_count_min, rem_ints, rem_exts;     ///< local variables and remainders
    int global_int_count, global_ext_count;                                  ///< global variables
    int global_nintci, global_nintcf, global_nextci, global_nextcf;

    // declare arrays for the complete arrays (global arrays to be distributed)
    int *lcc_linear, *elem_part, *epart;
    int **global_lcc;
    double *global_bs, *global_be, *global_bn, *global_bw, *global_bl, *global_bh, *global_bp, *global_su; 

    // declare arrays for data counts and displacements for each process in the workgroup
    int *local_int_count = calloc(nprocs, sizeof(int));
    int *local_ext_count = calloc(nprocs, sizeof(int));
    int *local_buff_int_count = malloc(nprocs * sizeof(int));
    int *global_int_start_idx = malloc(nprocs * sizeof(int));
    int *global_buff_int_start_idx = malloc(nprocs * sizeof(int));

    // MPI variables
    MPI_Request *requests;
    MPI_Status *statuses;

    // Initialize MPI requests and statuses depending on rank
    if(myrank == 0)
    {
        requests = (MPI_Request*)malloc((nprocs - 1) * sizeof(*requests));
        statuses = (MPI_Status*)malloc((nprocs - 1) * sizeof(*statuses));

        for(i = 0; i < (nprocs - 1); i++)
        {
            requests[i] = MPI_REQUEST_NULL;
        }
    }
    else
    {
        requests = (MPI_Request *)malloc(sizeof(*requests));
        statuses = (MPI_Status *)malloc(sizeof(*statuses));

        requests[0] = MPI_REQUEST_NULL;
    }

    // Allow only one process to read the file at a time
    if(myrank == 0)
    {
        //  read data for root only, nprocs set to 1
        f_status = read_binary_geo(file_name, nintci, nintcf, nextci, nextcf,
                                            &global_lcc,
                                            &global_bs, &global_be, &global_bn, &global_bw, 
                                            &global_bl, &global_bh, &global_bp, &global_su,
                                            points_count, points, elems);
        if ( f_status != 0 ) 
        {
            printf("ERROR: Data file not read.\n"); 
            return f_status;
        }
    }

    // broadcast all the indices to ALL procs from the root
    MPI_Bcast(nintci, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nintcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nextci, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(nextcf, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(points_count, 1, MPI_INT, 0, MPI_COMM_WORLD);

    global_int_count = (*nintcf) - (*nintci) + 1;
    global_ext_count = (*nextcf) - (*nextci) + 1;

    global_nintci = *nintci;
    global_nintcf = *nintcf;
    global_nextci = *nextci;
    global_nextcf = *nextcf;

    // Setup dummy elems and points for slave processes
    if(myrank != 0)
    {
        *elems = (int *)malloc(sizeof(int));
        *points = (int **)malloc((*points_count) * sizeof(int*));
        for(i = 0; i < *points_count; i++)
        {
            (*points)[i] = (int *)malloc(3 * sizeof(int));
        }
    }
    
    epart = (int *)malloc(global_int_count * sizeof(int));

    // Determine the type of data distribution method to implement
    if(!strcmp(part_type, "classic"))
    {
        if (myrank == 0)
        {
            lcc_linear = (int *) malloc((*nintcf + 1) * 6 * sizeof(int));
            linearizeArray(&global_lcc, lcc_linear, (*nintcf + 1));
        }

        // Determine int, double and points count for each process
        perProcDataDistribution(global_int_count, local_int_count, &local_int_count_min, &rem_ints, nprocs, myrank);
        
        // Determine subdomain internal starting indices
        for (i = 0; i < nprocs; i++)
        {
            global_int_start_idx[i] = i * local_int_count_min + (*nintci);
        }

        // Send the element partitions for index mapping
        elem_part = (int *)malloc(local_int_count[myrank] * sizeof(int));
        int **send_elem_part = (int **)malloc((nprocs - 1) * sizeof(int*));

        if(myrank == 0)
        {
            // First grab element partition for rank 0 (master rank)
            for(i = 0; i < local_int_count[myrank]; i++)
            {
                elem_part[i] = global_int_start_idx[myrank] + i;
            }

            // Send the remining ranks their element partitions
            for(i = 1; i < nprocs; i++)
            {
                // send_elem_part = (int *)malloc(local_int_count[i] * sizeof(int));
                send_elem_part[i - 1] = (int *)malloc(local_int_count[i] * sizeof(int));
                
                for(j = 0; j < local_int_count[i]; j++)
                {
                    // send_elem_part[j] = global_int_start_idx[i] + j;
                    send_elem_part[i - 1][j] = global_int_start_idx[i] + j;
                }

                MPI_Isend(send_elem_part[i - 1], local_int_count[i], MPI_INT, i, 0, MPI_COMM_WORLD, &(requests[i - 1]));
            }

            // Specify the element partition for each processor
            cnt = local_int_count_min;
            for(i = 0; i < nprocs; i++)
            {
                // Change the count for the last process to get remainder elements
                if(i == nprocs - 1)
                {
                    cnt = local_int_count[i];
                }

                for(j = 0; j < cnt; j++)
                {
                    epart[i*local_int_count_min + j] = i;
                }
            }
        }
        else
        {
            MPI_Irecv(elem_part, local_int_count[myrank], MPI_INT, 0, 0, MPI_COMM_WORLD, &(requests[0]));
        }

        // alloc memory for local lcc arrays
        int *lcc_linear_local;
        lcc_linear_local = (int *) calloc(local_int_count[myrank] * 6, sizeof(int));

        // set buffer sizes/counts for linearized arrays
        for (i = 0; i < nprocs; i++)
        {
            local_buff_int_count[i] = local_int_count[i] * 6;
            global_buff_int_start_idx[i] = global_int_start_idx[i] * 6;
        }

        // scatter lcc
        MPI_Scatterv(lcc_linear, local_buff_int_count, global_buff_int_start_idx, MPI_INT, lcc_linear_local, local_buff_int_count[myrank], MPI_INT, 0, MPI_COMM_WORLD);

        // allocate memory for each proc based on the above distribution
        *lcc = (int**) malloc(local_int_count[myrank] * sizeof(int*));
        for(i = 0; i < local_int_count[myrank]; i++)
        {
            (*lcc)[i] = (int *) malloc(6 * sizeof(int));
        }

        // re-structure linear lcc to lcc-structure of the code and get local double count
        unLinearizeArray(lcc, lcc_linear_local, local_int_count[myrank], global_int_count, &(local_ext_count[myrank]));

        // Send epart to all processes
        MPI_Bcast(epart, global_int_count, MPI_INT, 0, MPI_COMM_WORLD);

        find_neighbor_procs(*lcc, epart, nghb_cnt, nghb_to_rank, global_nintcf, nprocs, local_int_count[myrank], myrank);

        // Waits for the previous Isends (rank = 0) or Irecv (rank != 0) to complete
        if(myrank != 0)
        {
            MPI_Wait(&(requests[0]), &(statuses[0]));   
        }

        // Map local to global indices
        map_local_global_indices(local_global_index, elem_part, local_int_count[myrank]);

        // find number of cells to send to current proc's neighbours
        recvSendCountList(*lcc, epart, (*nghb_cnt), *nghb_to_rank, send_cnt, recv_cnt, send_lst, recv_lst, 
                    local_int_count[myrank], myrank, global_nintcf, *local_global_index, *global_local_index);

        // sort
        sortListForAll((*nghb_cnt), *recv_cnt, recv_lst, myrank, *global_local_index);
        sortListForAll((*nghb_cnt), *send_cnt, send_lst, myrank, *global_local_index);
        
        // Map global to local indices
        map_global_local_indices(global_local_index, epart, *lcc, *recv_cnt, *recv_lst, *nghb_cnt, global_int_count, global_ext_count, local_int_count[myrank], global_nextci, myrank);

        // Set local internal and external index ranges
        *nintci = 0;
        *nintcf = local_int_count[myrank] - 1;
        *nextci = local_int_count[myrank];
        *nextcf = local_ext_count[myrank] + local_int_count[myrank] - 1;

        // allocate local memory for B's
        *bs = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *be = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *bn = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *bw = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *bl = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *bh = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *bp = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *su = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));

        // copy doubles to all procs
        MPI_Scatterv(global_bs, local_int_count, global_int_start_idx, MPI_DOUBLE, bs[0], local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(global_be, local_int_count, global_int_start_idx, MPI_DOUBLE, be[0], local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(global_bn, local_int_count, global_int_start_idx, MPI_DOUBLE, bn[0], local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(global_bw, local_int_count, global_int_start_idx, MPI_DOUBLE, bw[0], local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(global_bl, local_int_count, global_int_start_idx, MPI_DOUBLE, bl[0], local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(global_bh, local_int_count, global_int_start_idx, MPI_DOUBLE, bh[0], local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(global_bp, local_int_count, global_int_start_idx, MPI_DOUBLE, bp[0], local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(global_su, local_int_count, global_int_start_idx, MPI_DOUBLE, su[0], local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        free(elem_part);
        free(lcc_linear_local);

        for(i = 0; i < (nprocs - 1); i++)
        {
            free(send_elem_part[i]);
        }
        free(send_elem_part);
    }
    else
    {
        double *bs_linear, *be_linear, *bn_linear, *bw_linear, *bl_linear, *bh_linear, *bp_linear,
                *su_linear;
        double *bs_linear_local, *be_linear_local, *bn_linear_local, *bw_linear_local, *bl_linear_local, *bh_linear_local, *bp_linear_local,
                *su_linear_local;
        
        idx_t *epart_metis = (idx_t *)malloc( (idx_t)global_int_count * sizeof(idx_t));

        // root reads all the data and calls the metis routines
        if (myrank == 0)
        {
            f_status = metis_partition(part_type, global_int_count, *points_count, nprocs, *elems, &epart_metis);

            if ( f_status != 0 ) return f_status;

            elem_part = (int *)malloc(global_int_count * nprocs * sizeof(int));

            for (i = 0; i < global_int_count; i++)
            {
                epart[i] = epart_metis[i];
            }

            for (j = 0; j < nprocs; j++)
            {
                local_int_count[j] = 0;
                for(i = 0; i < global_int_count; i++)
                {
                    if(epart[i] == j)
                    {
                        elem_part[ (j * global_int_count) + (local_int_count[j]++) ] = i;
                    }
                }
            }

            lcc_linear = (int *) malloc(global_int_count * 6 * sizeof(int));
            bs_linear = (double *) malloc((global_nextcf + 1) * sizeof(double));
            be_linear = (double *) malloc((global_nextcf + 1) * sizeof(double));
            bn_linear = (double *) malloc((global_nextcf + 1) * sizeof(double));
            bw_linear = (double *) malloc((global_nextcf + 1) * sizeof(double));
            bl_linear = (double *) malloc((global_nextcf + 1) * sizeof(double));
            bh_linear = (double *) malloc((global_nextcf + 1) * sizeof(double));
            bp_linear = (double *) malloc((global_nextcf + 1) * sizeof(double));
            su_linear = (double *) malloc((global_nextcf + 1) * sizeof(double));
        }

        // Send epart to all processes
        MPI_Bcast(epart, global_int_count, MPI_INT, 0, MPI_COMM_WORLD);

        // send local counts and displacements for each process
        MPI_Bcast(local_int_count, nprocs, MPI_INT, 0, MPI_COMM_WORLD);

        if (myrank == 0)
        {
            int local_int_end_index = 0;
            for (j = 0; j < nprocs; j++)
            {
                if (j == 0)
                {
                    global_int_start_idx[0] = (*nintci);                    ///< root start index is always zero or nintci
                    local_int_end_index = global_int_start_idx[j] + local_int_count[j];
                }
                else
                {
                    global_int_start_idx[j] = global_int_start_idx[j - 1] + local_int_count[j - 1];  //< Determine subdomain internal starting indices
                    local_int_end_index = global_int_start_idx[j] + local_int_count[j] + 1;
                    
                    // send element partition to all procs including root
                    MPI_Send(&elem_part[j * global_int_count], global_int_count, MPI_INT, j, 0, MPI_COMM_WORLD);
                }

                local_buff_int_count[j] = local_int_count[j] * 6;           ///< set buffer sizes/counts for linearized arrays
                global_buff_int_start_idx[j] = global_int_start_idx[j] * 6; ///< set buffer displacement for linearized arrays

                // Map local to global indices
                map_local_global_indices(local_global_index, &elem_part[j * global_int_count], local_int_count[j]);

                local_ext_count[j] = 0;

                int offset = global_int_start_idx[j];
                for (i = 0; i < local_int_count[j]; i++)
                {
                    // linearize global lcc using element index mapping
                    lcc_linear[( (i + offset) * 6) + 0] = global_lcc[ (*local_global_index)[i] ][0];
                    lcc_linear[( (i + offset) * 6) + 1] = global_lcc[ (*local_global_index)[i] ][1];
                    lcc_linear[( (i + offset) * 6) + 2] = global_lcc[ (*local_global_index)[i] ][2];
                    lcc_linear[( (i + offset) * 6) + 3] = global_lcc[ (*local_global_index)[i] ][3];
                    lcc_linear[( (i + offset) * 6) + 4] = global_lcc[ (*local_global_index)[i] ][4];
                    lcc_linear[( (i + offset) * 6) + 5] = global_lcc[ (*local_global_index)[i] ][5];
    
                    // Count the number of external cells
                    if(global_lcc[ (*local_global_index)[i] ][0] > global_nintcf) local_ext_count[j]++;
                    if(global_lcc[ (*local_global_index)[i] ][1] > global_nintcf) local_ext_count[j]++;
                    if(global_lcc[ (*local_global_index)[i] ][2] > global_nintcf) local_ext_count[j]++;
                    if(global_lcc[ (*local_global_index)[i] ][3] > global_nintcf) local_ext_count[j]++;
                    if(global_lcc[ (*local_global_index)[i] ][4] > global_nintcf) local_ext_count[j]++;
                    if(global_lcc[ (*local_global_index)[i] ][5] > global_nintcf) local_ext_count[j]++;

                    bs_linear[i + offset] = global_bs[ (*local_global_index)[i] ];
                    be_linear[i + offset] = global_be[ (*local_global_index)[i] ];
                    bn_linear[i + offset] = global_bn[ (*local_global_index)[i] ];
                    bw_linear[i + offset] = global_bw[ (*local_global_index)[i] ];
                    bl_linear[i + offset] = global_bl[ (*local_global_index)[i] ];
                    bh_linear[i + offset] = global_bh[ (*local_global_index)[i] ];
                    bp_linear[i + offset] = global_bp[ (*local_global_index)[i] ];
                    su_linear[i + offset] = global_su[ (*local_global_index)[i] ];
                }
        
                free(*local_global_index);
            } // end j-loop
        } // end root if
        else  // receive the element partitions
        {
            elem_part = (int *)malloc(global_int_count * sizeof(int));

            MPI_Irecv(elem_part, global_int_count, MPI_INT, 0, 0, MPI_COMM_WORLD, &(requests[0]));
        }

        // send local counts and displacements for each process
        MPI_Bcast(local_ext_count, nprocs, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(local_buff_int_count, nprocs, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(global_int_start_idx, nprocs, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(global_buff_int_start_idx, nprocs, MPI_INT, 0, MPI_COMM_WORLD);

        // allocate heap memory for the linearized local B's
        bs_linear_local = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        be_linear_local = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        bn_linear_local = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        bw_linear_local = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        bl_linear_local = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        bh_linear_local = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        bp_linear_local = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        su_linear_local = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));

        // copy doubles to all procs
        MPI_Scatterv(bs_linear, local_int_count, global_int_start_idx, MPI_DOUBLE, bs_linear_local, local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(be_linear, local_int_count, global_int_start_idx, MPI_DOUBLE, be_linear_local, local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(bn_linear, local_int_count, global_int_start_idx, MPI_DOUBLE, bn_linear_local, local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(bw_linear, local_int_count, global_int_start_idx, MPI_DOUBLE, bw_linear_local, local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(bl_linear, local_int_count, global_int_start_idx, MPI_DOUBLE, bl_linear_local, local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(bh_linear, local_int_count, global_int_start_idx, MPI_DOUBLE, bh_linear_local, local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(bp_linear, local_int_count, global_int_start_idx, MPI_DOUBLE, bp_linear_local, local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatterv(su_linear, local_int_count, global_int_start_idx, MPI_DOUBLE, su_linear_local, local_int_count[myrank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // alloc memory for local lcc arrays
        int *lcc_linear_local;
        lcc_linear_local = (int *) calloc(local_int_count[myrank] * 6, sizeof(int));

        MPI_Scatterv(lcc_linear, local_buff_int_count, global_buff_int_start_idx, MPI_INT, lcc_linear_local, local_buff_int_count[myrank], MPI_INT, 0, MPI_COMM_WORLD);
        
        // allocate memory for each proc based on the above distribution
        *lcc = (int**) malloc(local_int_count[myrank] * sizeof(int*));
        for(i = 0; i < local_int_count[myrank]; i++)
        {
            (*lcc)[i] = (int *) malloc(6 * sizeof(int));
        }

        *bs = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *be = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *bn = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *bw = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *bl = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *bh = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *bp = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        *su = (double *) malloc((local_int_count[myrank] + local_ext_count[myrank]) * sizeof(double));
        
        // Waits for the previous Isends (rank = 0) or Irecv (rank != 0) to complete
        if(myrank != 0)
        {
            MPI_Wait(&(requests[0]), &(statuses[0]));
        }

        // Map local to global indices
        map_local_global_indices(local_global_index, elem_part, local_int_count[myrank]);

        // re-structure lcc_linear to original form
        for(i = 0; i < local_int_count[myrank]; i++)
        {
            (*lcc)[ i ][0] = lcc_linear_local[(i * 6) + 0];
            (*lcc)[ i ][1] = lcc_linear_local[(i * 6) + 1];
            (*lcc)[ i ][2] = lcc_linear_local[(i * 6) + 2];
            (*lcc)[ i ][3] = lcc_linear_local[(i * 6) + 3];
            (*lcc)[ i ][4] = lcc_linear_local[(i * 6) + 4];
            (*lcc)[ i ][5] = lcc_linear_local[(i * 6) + 5];

            (*bs)[ i ] = bs_linear_local[i];
            (*be)[ i ] = be_linear_local[i];
            (*bn)[ i ] = bn_linear_local[i];
            (*bw)[ i ] = bw_linear_local[i];
            (*bl)[ i ] = bl_linear_local[i];
            (*bh)[ i ] = bh_linear_local[i];
            (*bp)[ i ] = bp_linear_local[i];
            (*su)[ i ] = su_linear_local[i];
        }

        find_neighbor_procs(*lcc, epart, nghb_cnt, nghb_to_rank, global_nintcf, nprocs, local_int_count[myrank], myrank);

        // Set local internal and external index ranges
        *nintci = 0;
        *nintcf = local_int_count[myrank] - 1;
        *nextci = local_int_count[myrank];
        *nextcf = local_ext_count[myrank] + local_int_count[myrank] - 1;

        (*local_int_count_out) = local_int_count[myrank];
        (*local_ext_count_out) = local_ext_count[myrank];
        
        // find number of cells to send to current proc's neighbours
        recvSendCountList(*lcc, epart, (*nghb_cnt), *nghb_to_rank, send_cnt, recv_cnt, send_lst, recv_lst, 
                    local_int_count[myrank], myrank, global_nintcf, *local_global_index, *global_local_index);

        // sort
        sortListForAll((*nghb_cnt), *recv_cnt, recv_lst, myrank, *global_local_index);
        sortListForAll((*nghb_cnt), *send_cnt, send_lst, myrank, *global_local_index);

        // Map global to local indices
        map_global_local_indices(global_local_index, epart, *lcc, *recv_cnt, *recv_lst, *nghb_cnt, global_int_count, global_ext_count, local_int_count[myrank], global_nextci, myrank);

        if (myrank == 0)
        {
            free(bs_linear);
            free(be_linear);
            free(bn_linear);
            free(bw_linear);
            free(bl_linear);
            free(bh_linear);
            free(bp_linear);
            free(su_linear);
        }

        free(lcc_linear_local);
        free(bs_linear_local);
        free(be_linear_local);
        free(bn_linear_local);
        free(bw_linear_local);
        free(bl_linear_local);
        free(bh_linear_local);
        free(bp_linear_local);
        free(su_linear_local);
        free(elem_part);

        free(epart_metis);
    }

    // only root de-allocs because read_go_bin called by root only
    if (myrank == 0)
    {
        free(global_bs);
        free(global_be);
        free(global_bn);
        free(global_bw);
        free(global_bl);
        free(global_bh);
        free(global_bp);
        free(global_su);
        
        // free lcc
        for(i = 0; i < (*nintcf) + 1; i++)
        {
            free(global_lcc[i]);
        }
        free(global_lcc);
        free(lcc_linear);
    }

    free(local_int_count);
    free(local_ext_count);
    free(local_buff_int_count);
    free(global_int_start_idx);
    free(global_buff_int_start_idx);
    free(epart);
    free(requests);
    free(statuses);
    
    return 0;
}

void find_neighbor_procs(int **lcc, int *epart, int *nghb_cnt, int **nghb_to_rank, int global_nintcf, int nprocs,
    int local_int_count, int myrank)
{
     int i, j, k, is_new_rank, elem_rank;

    // Initialize neighbor-to-rank array
    *nghb_to_rank = (int *)malloc((nprocs - 1) * sizeof(int));
    for (i = 0; i < (nprocs - 1); i++)
    {
        (*nghb_to_rank)[i] = -1;
    }

    // Check for local llc values located in neighboring processes
    *nghb_cnt = 0;
    for(i = 0; i < local_int_count; i++)
    {
        for (j = 0; j < 6; j++)
        {
            elem_rank = epart[lcc[i][j]];

            if((lcc[i][j] <= global_nintcf) && (elem_rank != myrank))
            {
                is_new_rank = TRUE;

                for(k = 0; k < (nprocs - 1); k++)
                {
                    if((*nghb_to_rank)[k] != -1)
                    {
                        if(elem_rank == (*nghb_to_rank)[k])
                        {
                            is_new_rank = FALSE;
                            break;
                        }
                    }
                    else
                    {
                        break;
                    }
                }

                if(is_new_rank == TRUE)
                {
                    (*nghb_to_rank)[(*nghb_cnt)++] = elem_rank;
                }
            }
        }
    }
}

int array_search(int *array, const int val, const int length)
{
    int i;
    int match = FALSE;

    for (i = 0; i < length; i++)
    {
        if ( array[i] == val )
        {
            match = TRUE;
            break;
        }
    }

    return match;
}

void recvSendCountList(int **lcc, int *epart, int nghb_cnt, int *nghb_to_rank, 
                    int **send_cnt, int **recv_cnt, int ***send_lst, int ***recv_lst, 
                    int local_int_count, int myrank, int global_nintcf, 
                    int *local_global_index, int *global_local_index)
{
    int i, j, k, neighbour, ret;

    (*send_cnt) = (int *) malloc(sizeof(int) * nghb_cnt);
    (*recv_cnt) = (int *) malloc(sizeof(int) * nghb_cnt);
    *send_lst = (int**) calloc(nghb_cnt, sizeof(int*));
    *recv_lst = (int**) calloc(nghb_cnt, sizeof(int*));

    // find send cell list for each neighbour by looping through nghb_cnt
    for (k = 0; k < nghb_cnt; k++)
    {
        (*send_cnt)[k] = 0;
        (*recv_cnt)[k] = 0;
        (*send_lst)[k] = (int*) malloc(local_int_count * 6 * sizeof(int));
        (*recv_lst)[k] = (int*) malloc(local_int_count * 6 * sizeof(int));

        for (i = 0; i < (local_int_count * 6); i++)
        {
            (*send_lst)[k][i] = -1;
            (*recv_lst)[k][i] = -1;
        }

        for (i = 0; i < local_int_count; i++)
        {
            for (j = 0; j < 6; j++)
            {
                neighbour = lcc[i][j];
                if ( (neighbour <= global_nintcf) )
                {
                    if ( (epart[ neighbour ] == nghb_to_rank[k]) )
                    {
                        ret = array_search( (*send_lst)[k], local_global_index[i], (*send_cnt)[k]);
                        if (ret == 0)
                        {
                            (*send_lst)[k][ (*send_cnt)[k] ] = local_global_index[i];
                            (*send_cnt)[k]++;
                        }

                        ret = array_search( (*recv_lst)[k], neighbour, (*recv_cnt)[k]);
                        if (ret == 0)
                        {
                            (*recv_lst)[k][ (*recv_cnt)[k] ] = neighbour;
                            (*recv_cnt)[k]++;
                        }
                    }
                }
            }
        }
    }
}

void sortListForAll(int nghb_cnt, int *sendRecv_cnt, int ***sendRecv_lst, 
                    int myrank, int *global_local_index)
{
    int i, k;

    for (k = 0; k < nghb_cnt; k++)
    {
        sortGlobalListGlobally((*sendRecv_lst)[k], sendRecv_cnt[k], myrank);
    }
}

void sortGlobalListGlobally(int *list, int count, int myrank)
{
    int i, local_val;
    int *temp_global_list;
    
    // the global list is sorted based on this local list
    temp_global_list = (int *) malloc(sizeof(int) * count);

    for (i = 0; i < count; i++)
    {
        temp_global_list[i] = list[i];
    }

    // sort the local list
    mergeSort(temp_global_list, count);
    
    for (i = 0; i < count; i++)
    {
        list[i] = temp_global_list[i];
    }

    free(temp_global_list);
}

void mergeSort(int *array, int length)
{
    int i, lengthLeft;
    int *arrayLeft, *arrayRight;

    if (length > 1)
    {
        // Divide array into 2 parts 
        lengthLeft = length/2;

        // allocate memory for both arrays. array = [arrayLeft + arrayRight]
        arrayLeft   = (int *) malloc(sizeof(int) * lengthLeft);
        arrayRight  = (int *) malloc(sizeof(int) * (length - lengthLeft));

        // copy array into left and right parts
        for (i = 0; i < lengthLeft; i++)
        {
            arrayLeft[i] = array[i];
        }

        for (i = 0; i < (length - lengthLeft); i++)
        {
            arrayRight[i] = array[lengthLeft + i];
        }

        // recursive call
        mergeSort(arrayLeft, lengthLeft);
        mergeSort(arrayRight, (length - lengthLeft));

        // merge the sorted arrays
        mergeSortedArrays(arrayLeft, arrayRight, array, lengthLeft, (length - lengthLeft), length);

        free(arrayLeft);
        free(arrayRight);
    }
}

void mergeSortedArrays(int *arrayLeft, int *arrayRight, int *array, int lengthLeft, int lengthRight, int length)
{
    int k;
    int i = 0;
    int j = 0;

    for (k = 0; k < length; k++)
    {
        // if left array complete then add remaining right elems to array
        if ( i > (lengthLeft - 1) )
        {
            array[k] = arrayRight[j];
            j++;
        }
        else if ( j > (lengthRight - 1) )
        {
            array[k] = arrayLeft[i];
            i++;
        }
        else if ( arrayLeft[i] < arrayRight[j] )
        {
            array[k] = arrayLeft[i];
            i++;
        }
        else
        {
            array[k] = arrayRight[j];
            j++;
        }
    }
}

int metis_partition(char *part_type, int global_int_count, int points_count, int nprocs, int *elems, idx_t **epart)
{
    int f_status, i;

    // METIS parameters
    idx_t *eptr, *eind, *objval, *npart;
    idx_t ne = (idx_t)global_int_count;
    idx_t nn = (idx_t)points_count;
    idx_t nparts = (idx_t)nprocs;

    // Allocate memory and initialize METIS structures
    eptr = (idx_t *)malloc((global_int_count + 1) * sizeof(idx_t));
    eind = (idx_t *)malloc(global_int_count * 8 * sizeof(idx_t));
    npart = (idx_t *)malloc(points_count * sizeof(idx_t));
    objval = (idx_t *)malloc( sizeof(idx_t) );

    for(i = 0; i < (global_int_count + 1); i++)
    {
        eptr[i] = (idx_t)(8 * i);
    }

    for(i = 0; i < (global_int_count * 8); i++)
    {
        eind[i] = (idx_t)elems[i];
    }

    if(!strcmp(part_type, "dual"))
    {
        idx_t ncommon = (idx_t)4;
        f_status = METIS_PartMeshDual(  &ne, &nn, eptr, eind, NULL, NULL, &ncommon, &nparts, NULL, NULL,
                                         objval, *epart, npart);
    }
    else if(!strcmp(part_type, "nodal"))
    {
        f_status = METIS_PartMeshNodal( &ne, &nn, eptr, eind, NULL, NULL, &nparts, NULL, NULL,
                                         objval, *epart, npart);
    }
    else
    {
        printf("ERROR: The METIS partition type %s is not supported\n", part_type);
        return 1;
    }
    if(f_status != METIS_OK)
    {
        printf("ERROR: METIS returned something other than METIS_OK\n");
        return 1;
    }
    
    // Deallocate memory
    free(eptr);
    free(eind);
    free(npart);
    free(objval);

    return 0;
}

void linearizeArray(int ***lcc_bad, int *lcc_good, const int local_int_count)
{
    int i, j;

    for (i = 0; i < local_int_count; i++)
    {
        for (j = 0; j < 6; j++)
        {
            lcc_good[(i * 6) + j] = (*lcc_bad)[i][j];
        }
    }
}

void unLinearizeArray(int ***lcc_bad, int *lcc_good, const int local_int_count, const int global_int_count, int *local_ext_count)
{
    int i, j;

    for (i = 0; i < local_int_count; i++)
    {
        for (j = 0; j < 6; j++)
        {
            (*lcc_bad)[i][j] = lcc_good[(i * 6) + j];

            // Count the number of external (double) cells
            if((*lcc_bad)[i][j] > global_int_count) (*local_ext_count)++;
        }
    }
}

void perProcDataDistribution(int totalCount, int *perProcCount, int *quotient, int *rem, int numProcs, int myrank)
{
    (*quotient) = totalCount / numProcs;        ///< what each proc will definitely get
    (*rem) = totalCount % numProcs;             ///< remainder is added to this, if needed

    int i;
    for (i = 0; i < numProcs; i++)
    {
        perProcCount[i] = (*quotient);
    }

    perProcCount[numProcs - 1] += (*rem);
}

void map_local_global_indices(int **local_global_index, int *elem_part, int part_count)
{
    int i;

    // Allocate memory
    *local_global_index = malloc(part_count * sizeof(int));

    // Construct mapping
    for(i = 0; i < part_count; i++)
    {
        (*local_global_index)[i] = elem_part[i];
    }
}

void map_global_local_indices(int **global_local_index, int *epart, int **lcc, int *recv_cnt, int **recv_lst,
    int nghb_cnt, int global_int_count, int global_ext_count, int local_int_count, int global_nextci, int myrank)
{
    int i, j;
    int k = 0;
    int global_count = global_int_count + global_ext_count;

    // Allocate memory
    *global_local_index = (int *)malloc(global_count * sizeof(int));

    // Initialize mapping elements
    for(i = 0; i < global_count; i++)
    {
        (*global_local_index)[i] = -1;
    }

    // Map internal cells
    for(i = 0; i < global_int_count; i++)
    {
        if(epart[i] == myrank)
        {
            (*global_local_index)[i] = k++;
        }
    }

    // Map external cells
    for(i = 0; i < local_int_count; i++)
    {
        for(j = 0; j < 6; j++)
        {
            if(lcc[i][j] >= global_nextci)
            {
                (*global_local_index)[lcc[i][j]] = k++;
            }
        }
    }

    // Map ghost cells
    for(i = 0; i < nghb_cnt; i++)
    {
        for(j = 0; j < recv_cnt[i]; j++)
        {
            (*global_local_index)[recv_lst[i][j]] = k++;
        }
    }
}
