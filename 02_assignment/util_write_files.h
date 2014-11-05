/**
 * Helper functions for writing results to VTK and text files
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#ifndef XWRITE_H_
#define XWRITE_H_

/**
 * Export the simulation statistics to a text file
 *
 * @param in_file_name input file used for this simulation
 * @param out_file_name target output file
 * @param nintci index of the first element to include in the output file
 * @param nintcf index of the last element to include in the output file
 * @param var array with the computed results
 * @param total_iters total number of iterations required to compute the solution
 * @param residual_ratio the final residual ratio
 * @return
 */
int store_simulation_stats(char *in_file_name, char *out_file_name, int nintci, int nintcf,
                           double *var, int total_iters, double residual_ratio);

/**
 * Write a VTK header with the desired geometry data to a file
 *
 * @param experiment_name title of the VTK file (in out case simply the used input file)
 * @param out_file_name target output file
 * @param start_index index of the first element to include in the output file
 * @param end_index index of the last element to include in the output file
 * @param points_count total number of points
 * @param points coordinates of the points that define the volume cells
 * @param elems neighboring information for the cells
 */
void vtk_write_unstr_grid_header(const char *experiment_name, const char *out_file_name,
                                 int start_index, int end_index, int points_count, int **points,
                                 int *elems);

/**
 * Append the values of a double variable to a VTK file
 *
 * @param out_file_name target output file
 * @param var_name name of the variable
 * @param start_index index of the first element to include in the output file
 * @param end_index index of the last element to include in the output file
 * @param values array with the values of var_name at the grid points
 */
void vtk_append_double(const char *out_file_name, const char *var_name, int start_index,
                       int end_index, double *values);

/**
 * Append the values of an integer variable to a VTK file
 *
 * @param out_file_name target output file
 * @param var_name name of the variable
 * @param start_index index of the first element to include in the output file
 * @param end_index index of the last element to include in the output file
 * @param values array with the values of var_name at the grid points
 */
void vtk_append_integer(const char *out_file_name, const char *var_name, int start_index,
                        int end_index, int *values);

#endif /* XWRITE_H_ */

