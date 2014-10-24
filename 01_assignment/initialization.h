/**
 * Initialization step - parse the input file, compute data distribution, initialize LOCAL computational arrays
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#ifndef INITIALIZATION_H_
#define INITIALIZATION_H_

int initialization(char* file_in, char* file_format, int* nintci, int* nintcf, int* nextci,
                   int* nextcf, int*** lcc, double** bs, double** be, double** bn, double** bw,
                   double** bl, double** bh, double** bp, double** su, double** var, double** cgup, 
                   double** oc, double** cnorm);

#endif /* INITIALIZATION_H_ */

