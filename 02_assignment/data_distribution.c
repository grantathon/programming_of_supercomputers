#include <stdlib.h>

#include "data_distribution.h"

int read_data(     char *file_name, char *part_type, char *read_type, int *nintci, int *nintcf, int *nextci, int *nextcf,
                                    int ***lcc, double **bs, double **be, double **bn, double **bw, double **bl,
                                    double **bh,double **bp, double **su, int *points_count, int ***points, int **elems,
                                    int **local_global_index, int nprocs, int myrank)
{
    int f_status = 0;

    if(!strcmp(read_type, "oneread"))
    {
        // f_status = one_proc_read();
    }
    else if(!strcmp(read_type, "allread"))
    {
        // f_status = all_proc_read();
    }
    else
    {
        if(myrank == 0)
        {
            printf("ERROR: Unsupported input file reading type: %s\n", read_type);  
        }

        return 1;
    }

    if(f_status == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int all_proc_read()
{

}

int one_proc_read()
{
    
}

int classic_partition(  char *file_name, char *read_type, int *nintci, int *nintcf, int *nextci, int *nextcf,
                                        int ***lcc, double **bs, double **be, double **bn, double **bw, double **bl,
                                        double **bh,double **bp, double **su, int *points_count, int ***points, int **elems,
                                        int **local_global_index, int nprocs, int myrank)
{
    int f_status = 0;

    // Read the input data based on the specified read type
    // f_status = read_data( file_name, read_type, nintci, nintcf, nextci, nextcf, lcc, bs, be, bn, bw, bl, bh,
    //                                     bp, su, points_count, points, elems, local_global_index, nprocs, myrank);

    if ( f_status != 0 ) return f_status;

    // Determine number of elements and ghosts to distribute
    int global_elem_count = (*nintcf) - (*nintci) + 1;
    int local_elem_count = global_elem_count / nprocs;
    int rem_elems = global_elem_count % nprocs;

    int global_ghost_count = (*nextcf) - (*nextci) + 1;
    int local_ghost_count = global_ghost_count / nprocs;
    int rem_ghosts = global_ghost_count % nprocs;

    // Determine subdomain (i.e., element and ghost indices) to assign to current processor
    *nintci = (*myrank)*local_elem_count;
    *nintcf = (*nintci) + local_elem_count - 1;

    *nextci = (*myrank)*local_ghost_count;
    *nextcf = (*nextci) + local_ghost_count - 1;

    // Add remaining elements and ghosts to last processes (if necessary)
    if(myrank == nprocs - 1 && rem_elems != 0)
    {
        local_elem_count += rem_elems;
        *nintcf += rem_elems;

        local_ghost_count += rem_ghosts;
        *nextcf += rem_ghosts;
    }

    return 0;
}

int partition_data(  char *file_name, char *part_type, int *nintci, int *nintcf, int *nextci, int *nextcf,
                                int ***lcc, double **bs, double **be, double **bn, double **bw, double **bl,
                                double **bh,double **bp, double **su, int *points_count, int ***points, int **elems,
                                int **local_global_index, int nprocs, int myrank)
{
    int f_status = 0;

    // Determine the type of data distribution method to implement
    if(!strcmp(part_type, "classic"))
    {
        f_status =  classic_partition(  file_name, read_type, nintci, nintcf, nextci, nextcf, lcc, bs, be, bn, bw, bl,
                                                         bh, bp, su, points_count, points, elems, local_global_index, nprocs, myrank);
    }
    else if(!strcmp(part_type, "dual"))
    {
        // f_status =  dual_partition(  file_name, read_type, nintci, nintcf, nextci, nextcf, lcc, bs, be, bn, bw, bl,
        //                                                  bh, bp, su, points_count, points, elems, local_global_index, nprocs, myrank);
    }
    else if(!strcmp(part_type, "nodal"))
    {
        // f_status =  nodal_partition( file_name, read_type, nintci, nintcf, nextci, nextcf, lcc, bs, be, bn, bw, bl,
        //                                                   bh, bp, su, points_count, points, elems, local_global_index, nprocs, myrank);
    }
    else
    {
        if(myrank == 0)
        {
            printf("ERROR: Unsupported data partition type: %s\n", part_type);  
        }

        return 1;
    }

    // f_status = map_local_global_indices();

    if(f_status == 0)
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

int map_local_global_indices()
{
    return 0;
}