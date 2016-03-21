/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012, 13-Nov-2014
 * @author V. Petkov, A. Berariu
 */

#ifndef INITIALIZATION_H_
#define INITIALIZATION_H_

int initialization(char* file_in, char* part_type, char* read_type, int nprocs, int myrank,
                   int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, int* points_count,
                   int*** points, int** elems, double** var, double** cgup, double** oc,
                   double** cnorm, int** local_global_index, int** global_local_index,
                   int *nghb_cnt, int** nghb_to_rank, int** send_cnt, int*** send_lst, 
                   int **recv_cnt, int*** recv_lst);

#endif /* INITIALIZATION_H_ */

