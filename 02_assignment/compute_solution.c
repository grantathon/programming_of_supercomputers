/**
 * Computational loop
 *
 * @file compute_solution.c
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

int compute_solution_parallel(int nprocs, int myrank, const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                              double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                              double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                              int* local_global_index, int* global_local_index, int nghb_cnt, 
                              int* nghb_to_rank, int* send_cnt, int** send_lst, int *recv_cnt, int** recv_lst)
{
    /** parameters used in gccg */
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int i, j;
    int nomax = 3;

    // Parameters used in computation loop
    double oc1 = 0.0;
    double oc2 = 0.0;
    double occ = 0.0;

    // Reference residuals
    double resref = 0.0;

    // MPI communication parameters
    int ghost_offset = 0;
    int local_ghost_count = 0;
    int **block_send_length, **block_send_stride, **block_recv_length, **block_recv_stride;
    MPI_Datatype *mpi_send_type;
    MPI_Datatype *mpi_recv_type;
    MPI_Status *statuses = (MPI_Status *)malloc(2 * nghb_cnt * sizeof(MPI_Status));
    MPI_Request *requests = (MPI_Request *)malloc(2 * nghb_cnt * sizeof(MPI_Request));

    for(nc = 0; nc < nghb_cnt; nc++)
    {
        local_ghost_count += recv_cnt[nc];
    }

    for(nc = 0; nc < (2 * nghb_cnt); nc++)
    {
        requests[nc] = MPI_REQUEST_NULL;
    }

    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (nintcf + 1));

    // initialize the reference residual
    for ( nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }

    // Sum residuals from all procs and distribute value
    MPI_Allreduce(MPI_IN_PLACE, &resref, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    resref = sqrt(resref);

    if ( resref < 1.0e-15 ) {
        fprintf(stderr, "Residue sum less than 1.e-15 - %lf\n", resref);
        return 0;
    }

    /** the computation vectors */
    double *direc1 = (double *) calloc(sizeof(double), (nextcf + 1 + local_ghost_count));
    double *direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *adxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *adxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor2 = (double *) calloc(sizeof(double), (nintcf + 1));

    // Allocate memory for construction of custom indexed types
    block_send_length = (int **)malloc(nghb_cnt * sizeof(int*));
    block_send_stride = (int **)malloc(nghb_cnt * sizeof(int*));
    block_recv_length = (int **)malloc(nghb_cnt * sizeof(int*));
    block_recv_stride = (int **)malloc(nghb_cnt * sizeof(int*));

    mpi_send_type = (MPI_Datatype *)malloc(nghb_cnt * sizeof(MPI_Datatype));
    mpi_recv_type = (MPI_Datatype *)malloc(nghb_cnt * sizeof(MPI_Datatype));

    // Create the MPI indexed data types
    for(i = 0; i < nghb_cnt; i++)
    {
        // Allocate memory based on number of blocks
        block_send_length[i] = (int *)malloc(send_cnt[i] * sizeof(int));
        block_send_stride[i] = (int *)malloc(send_cnt[i] * sizeof(int));

        // Set block lengths to one
        for(j = 0; j < send_cnt[i]; j++)
        {
            block_send_length[i][j] = 1;
        }

        for(j = 0; j < send_cnt[i]; j++)
        {
            block_send_stride[i][j] = global_local_index[send_lst[i][j]];
        }

        // Allocate memory for new MPI indexed type
        MPI_Type_indexed(send_cnt[i], block_send_length[i], block_send_stride[i], MPI_DOUBLE, &(mpi_send_type[i]));
        MPI_Type_commit(&(mpi_send_type[i]));

        block_recv_length[i] = (int *)malloc(sizeof(int));
        block_recv_stride[i] = (int *)malloc(sizeof(int));

        block_recv_length[i][0] = recv_cnt[i];
        block_recv_stride[i][0] = nextcf + 1 + ghost_offset;

        // Allocate memory for new MPI indexed send type
        MPI_Type_indexed(1, block_recv_length[i], block_recv_stride[i], MPI_DOUBLE, &(mpi_recv_type[i]));
        MPI_Type_commit(&(mpi_recv_type[i]));

        ghost_offset += recv_cnt[i];
    }


    while ( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/

        // update the old values of direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

        // Communicate ghost cells between neighbor procs
        communicate_ghost_cells(direc1, nghb_cnt, nghb_to_rank, recv_cnt, send_cnt, recv_lst, send_lst, myrank, global_local_index, mpi_send_type,
                                mpi_recv_type, requests, statuses);

        // compute new guess (approximation) for direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc2[nc] =   bp[nc] * direc1[nc]
                         - bs[nc] * direc1[global_local_index[lcc[nc][0]]]
                         - be[nc] * direc1[global_local_index[lcc[nc][1]]]
                         - bn[nc] * direc1[global_local_index[lcc[nc][2]]]
                         - bw[nc] * direc1[global_local_index[lcc[nc][3]]]
                         - bl[nc] * direc1[global_local_index[lcc[nc][4]]]
                         - bh[nc] * direc1[global_local_index[lcc[nc][5]]];
        }

        /********** END COMP PHASE 1 **********/

        /********** START COMP PHASE 2 **********/

        // execute normalization steps
        if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;

            for ( nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + direc2[nc] * adxor1[nc];
            }

            MPI_Allreduce(MPI_IN_PLACE, &occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

            oc1 = occ / cnorm[1];
            for ( nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
            }

            if1++;
        } else {
            if ( nor1 == 2 ) {
                oc1 = 0;
                occ = 0;

                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + direc2[nc] * adxor1[nc];
                }

                MPI_Allreduce(MPI_IN_PLACE, &occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                oc1 = occ / cnorm[1];
                oc2 = 0;
                occ = 0;
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + direc2[nc] * adxor2[nc];
                }

                MPI_Allreduce(MPI_IN_PLACE, &occ, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

                oc2 = occ / cnorm[2];
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                    direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                }

                if2++;
            }
        }

        // compute the new residual
        cnorm[nor] = 0;
        double omega = 0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }

        // Sum omegas from all procs and distribute value
        MPI_Allreduce(MPI_IN_PLACE, &(cnorm[nor]), 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, &omega, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        omega = omega / cnorm[nor];

        double res_updated = 0.0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            res_updated = res_updated + resvec[nc] * resvec[nc];
            var[nc] = var[nc] + omega * direc1[nc];
        }

        // Sum residuals from all procs and distribute value
        MPI_Allreduce(MPI_IN_PLACE, &res_updated, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        res_updated = sqrt(res_updated);

        *residual_ratio = res_updated / resref;

        // exit on no improvements of residual
        if ( *residual_ratio <= 1.0e-10 ) break;

        iter++;

        // prepare additional arrays for the next iteration step
        if ( nor == nomax ) {
            nor = 1;
        } else {
            if ( nor == 1 ) {
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if ( nor == 2 ) {
                    for ( nc = nintci; nc <= nintcf; nc++ ) {
                        dxor2[nc] = direc1[nc];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        
        /********** END COMP PHASE 2 **********/
    }

    // Deallocate memory
    free(direc1);
    free(direc2);
    free(adxor1);
    free(adxor2);
    free(dxor1);
    free(dxor2);
    free(resvec);

    for(i = 0; i < nghb_cnt; i++ )
    {
        free(block_send_length[i]);
        free(block_recv_length[i]);
    }
    free(block_send_length);
    free(block_recv_length);

    for(i = 0; i < nghb_cnt; i++ )
    {
        free(block_send_stride[i]);
        free(block_recv_stride[i]);
    }
    free(block_send_stride);
    free(block_recv_stride);

    free(mpi_send_type);
    free(mpi_recv_type);

    return iter;
}

void communicate_ghost_cells(   double *direc1, int nghb_cnt, int *nghb_to_rank, int *recv_cnt, int *send_cnt, int **recv_lst, 
                                int **send_lst, int myrank, int *global_local_index, MPI_Datatype *mpi_send_type,
                                MPI_Datatype *mpi_recv_type, MPI_Request *requests, MPI_Status *statuses)
{
    int i, j, local_idx;

    // send data from direc1 using MPI_data_type
    for (i = 0; i < nghb_cnt; i++)
    {
        MPI_Isend(direc1, 1, mpi_send_type[i], nghb_to_rank[i], 0, MPI_COMM_WORLD, &requests[i]);
    }

    // send data for direc1 using MPI_data_type
    for (i = 0; i < nghb_cnt; i++)
    {
        MPI_Irecv(direc1, 1, mpi_recv_type[i], nghb_to_rank[i], 0, MPI_COMM_WORLD, &requests[nghb_cnt + i]);
    }

    MPI_Waitall((2 * nghb_cnt), requests, statuses);
}

int compute_solution_serial(int nprocs, int myrank, const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                            double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                            double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                            int* local_global_index, int* global_local_index, int nghb_cnt, 
                            int* nghb_to_rank, int* send_cnt, int** send_lst, int *recv_cnt, int** recv_lst){
    /** parameters used in gccg */
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int nomax = 3;

    /** the reference residual */
    double resref = 0.0;

    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (nintcf + 1));

    // initialize the reference residual
    for ( nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }

    resref = sqrt(resref);
    if ( resref < 1.0e-15 ) {
        fprintf(stderr, "Residue sum less than 1.e-15 - %lf\n", resref);
        return 0;
    }

    /** the computation vectors */
    double *direc1 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *adxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *adxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor2 = (double *) calloc(sizeof(double), (nintcf + 1));

    while ( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/
        // update the old values of direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

        // compute new guess (approximation) for direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[nc][0]]
                         - be[nc] * direc1[lcc[nc][1]] - bn[nc] * direc1[lcc[nc][2]]
                         - bw[nc] * direc1[lcc[nc][3]] - bl[nc] * direc1[lcc[nc][4]]
                         - bh[nc] * direc1[lcc[nc][5]];
        }
        /********** END COMP PHASE 1 **********/

        /********** START COMP PHASE 2 **********/
        // execute normalization steps
        double oc1, oc2, occ;
        if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;

            for ( nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + direc2[nc] * adxor1[nc];
            }

            oc1 = occ / cnorm[1];
            for ( nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
            }

            if1++;
        } else {
            if ( nor1 == 2 ) {
                oc1 = 0;
                occ = 0;

                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + direc2[nc] * adxor1[nc];
                }

                oc1 = occ / cnorm[1];
                oc2 = 0;
                occ = 0;
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + direc2[nc] * adxor2[nc];
                }

                oc2 = occ / cnorm[2];
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                    direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                }

                if2++;
            }
        }

        // compute the new residual
        cnorm[nor] = 0;
        double omega = 0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }

        omega = omega / cnorm[nor];
        double res_updated = 0.0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            res_updated = res_updated + resvec[nc] * resvec[nc];
            var[nc] = var[nc] + omega * direc1[nc];
        }

        res_updated = sqrt(res_updated);
        *residual_ratio = res_updated / resref;

        // exit on no improvements of residual
        if ( *residual_ratio <= 1.0e-10 ) break;

        iter++;

        // prepare additional arrays for the next iteration step
        if ( nor == nomax ) {
            nor = 1;
        } else {
            if ( nor == 1 ) {
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if ( nor == 2 ) {
                    for ( nc = nintci; nc <= nintcf; nc++ ) {
                        dxor2[nc] = direc1[nc];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        /********** END COMP PHASE 2 **********/
    }

    free(direc1);
    free(direc2);
    free(adxor1);
    free(adxor2);
    free(dxor1);
    free(dxor2);
    free(resvec);

    return iter;
}
