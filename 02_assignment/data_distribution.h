#ifndef DATA_DISTRIBUTION_H_
#define DATA_DISTRIBUTION_H_

int read_data(  char *file_name, char *part_type, char *read_type, int *nintci, int *nintcf, int *nextci, int *nextcf,
                         int ***lcc, double **bs, double **be, double **bn, double **bw, double **bl,
                         double **bh,double **bp, double **su, int *points_count, int ***points, int **elems,
                         int **local_global_index, int nprocs, int myrank);

#endif /* DATA_DISTRIBUTION_H_ */
