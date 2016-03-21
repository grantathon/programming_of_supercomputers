/*
 * compute_solution.h
 *
 *  Created on: Oct 21, 2012
 *      Author: petkovve
 */

#ifndef COMPUTE_SOLUTION_H_
#define COMPUTE_SOLUTION_H_

int compute_solution_parallel(int nprocs, int myrank, const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                              double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                              double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                              int* local_global_index, int* global_local_index, int nghb_cnt, 
                              int* nghb_to_rank, int* send_cnt, int** send_lst, int *recv_cnt, int** recv_lst);

int compute_solution_serial(int nprocs, int myrank, const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                            double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                            double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                            int* local_global_index, int* global_local_index, int nghb_cnt, 
                            int* nghb_to_rank, int* send_cnt, int** send_lst, int *recv_cnt, int** recv_lst);

#endif /* COMPUTE_SOLUTION_H_ */

