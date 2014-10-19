/**
 * Helper functions for writing results to VTK and text files
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#ifndef XWRITE_H_
#define XWRITE_H_

int write_result(char *infilename, char *outfilename, int nintci, int nintcf,
        double *var, int iter, double ratio);

int write_result_vtk(char *outFileName, int startIdx, int endIdx, int nodeCnt, 
	int **points, int **elems, double *vector);


#endif /* XWRITE_H_ */

