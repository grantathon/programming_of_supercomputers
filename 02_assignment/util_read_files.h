/**
 * Helper functions for reading from input data file
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#ifndef XREAD_H_
#define XREAD_H_

int read_binary_geo(char *file_name, int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF, int ***LCC,
                    double **BS, double **BE, double **BN, double **BW, double **BL, double **BH,
                    double **BP, double **SU, int* nodeCnt, int***points, int** elems);

#endif /* XREAD_H_ */


