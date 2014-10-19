/**
 * Finalization step - write results and other computational vectors to files
 *
 * @date 22-Oct-2012
 * @author V. Petkov
 */

#include <stdio.h>
#include "util_write_files.h"

void finalization(char* file_in, int total_iters, double residual_ratio,
                  int nintci, int nintcf, double* var, double* cgup, double* su){

    char file_out[30] = "summary.out";
    
    int status = write_result(file_in, file_out, nintci, nintcf, var, total_iters, residual_ratio);
  
    if ( status != 0 ) fprintf(stderr, "Error when trying to write to file %s\n", file_out);
}

