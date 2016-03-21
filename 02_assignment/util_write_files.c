/**
 * Helper functions for writing results to VTK and text files
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "util_write_files.h"
#include "test_functions.h"

int store_simulation_stats(char *in_file_name, char *out_file_name, int nintci, int nintcf,
                           double *var, int total_iters, double residual_ratio)
{
    double *points = (double *) malloc((nintcf + 1) * sizeof(double));
    int i1, i2, i3, i4, i5;

    int counter;
    for ( counter = nintci; counter <= nintcf; counter++ )
        points[counter] = counter;

    if ( nintcf <= 1 ) {
        fprintf(stderr, "Error: NINTCF <= 1\n");
        return -1;
    }

    i1 = nintcf + 1;

    while ( i1 != 0 ) {
        i1 = i1 / 2;
        i2 = nintcf + 1 - i1;
        i4 = 1;

        do {
            i3 = i4;
            do {
                i5 = i3 + i1;
                if ( var[i3] <= var[i5] ) break;

                double z_dum = var[i3], i_dum = points[i3];

                var[i3] = var[i5];
                points[i3] = points[i5];
                var[i5] = z_dum;
                points[i5] = i_dum;
                i3 = i3 - i1;
            } while ( i3 >= 1 );
            i4++;
        } while ( i4 < i2 );
    }

    FILE *fp = fopen(out_file_name, "w");
    if ( fp == NULL ) {
        fprintf(stderr, "Error opening file %s for writing\n", out_file_name);
        return -1;
    }

    printf("========================================\n");
    printf("= AVL -  Linear Equation Solver - GCCG =\n");
    printf("========================================\n\n");
    printf("Input File:  %s\n", in_file_name);
    printf("Output File:  %s\n", out_file_name);
    printf("No. of Active Cells:  %d\n", nintcf + 1);
    printf("Iterations Count: %d\n", total_iters);
    printf("Residual Ratio: %e\n", residual_ratio);
    printf("========================================\n\n");

    fprintf(fp, "========================================\n");
    fprintf(fp, "= AVL -  Linear Equation Solver - GCCG =\n");
    fprintf(fp, "========================================\n");
    fprintf(fp, "\n\n");
    fprintf(fp, "Input File:  %s\n", in_file_name);
    fprintf(fp, "========================================\n\n");
    fprintf(fp, "Output File:  %s\n", out_file_name);
    fprintf(fp, "========================================\n\n");
    fprintf(fp, "No. of Active Cells:  %d\n", nintcf + 1);
    fprintf(fp, "========================================\n\n");
    fprintf(fp, "Iterations Count: %d\n", total_iters);
    fprintf(fp, "========================================\n\n");
    fprintf(fp, "Residual Ratio: %e\n", residual_ratio);
    fprintf(fp, "========================================\n\n");
    fprintf(fp, "Addresses Solution (Minima)  \t Addresses Solution (Maxima)\n");
    fprintf(fp, "===========================         ===========================\n");

    int N;
    for ( N = 1; N <= 10; N++ )
        fprintf(fp, "%8.0lf \t %lf \t\t %8.0lf \t %lf\n", points[N], var[N], points[nintcf - N + 1],
                var[nintcf - N + 1]);

    fprintf(fp, "========================================\n");

    fclose(fp);

    free(points);
    return 0;
}

void vtk_write_unstr_grid_header(const char *experiment_name, const char *out_file_name,
                                 int start_index, int end_index, int points_count, int **points,
                                 int *elems) {
    int i, j;
    FILE *fp = NULL;
    // Total number of elements (cells)
    int elem_count = end_index - start_index + 1;

    fp = fopen(out_file_name, "w");
    if ( fp == NULL ) {
        fprintf(stderr, "Failed to open %s", out_file_name);
        return;
    }

    fprintf(fp, "# vtk DataFile Version 3.1\n");
    fprintf(fp, "%s\n", experiment_name);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    /*
     * The first line notes how many points there will be and in which format they'll be supplied (here "DOUBLE").
     * Each following line contains the xyz-coordinates of a point.
     * Based on this list, each point is assigned an id. The first point has id 0, the second point id 1, and so forth.
     */
    fprintf(fp, "POINTS %d DOUBLE\n", points_count);
    for ( i = 0; i < points_count; i++ )
        fprintf(fp, "%d %d %d\n", points[i][0], points[i][1], points[i][2]);

    fprintf(fp, "\n");

    /*
     * The first line notes how many cells there will be and how many numbers total will be supplied in the CELLS-block.
     * Each following cell line starts with a number saying how many point IDs are to be read in that line (here 8) followed by the list of those point IDs.
     * Based on this list, each cell is assigned an id. The first cell has id 0, the second id 1, etc.
     */
    fprintf(fp, "CELLS %d %d\n", elem_count, elem_count * 9);
    for ( i = 0; i < elem_count; i++ ) {
        fprintf(fp, "8 ");
        for ( j = 0; j < 8; j++ )
            fprintf(fp, "%d ", elems[8 * i + j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    /*
     * The first line gives how many cell-types are to be set (= the number of cells given in CELLS).
     * The following line is a list of cell-types that are assigned to all cells. The first cell is of type "11", so is the second cell and so forth.
     * Cell type "11" is VTK_VOXEL (volume cell) and requires 8 point coordinates to define the  of each cell (given in CELLS).
     */
    fprintf(fp, "CELL_TYPES %d\n", elem_count);
    for ( i = 0; i < elem_count; i++ )
        fprintf(fp, "11 ");
    fprintf(fp, "\n\n");

    /*
     * Denotes the beginning of the datasets that will be assigned to the cell.
     * POINT_DATA can be used to assign the datasets to the points instead of the cells.
     */
    fprintf(fp, "CELL_DATA %d\n", elem_count);
    fprintf(fp, "\n");

    if ( fclose(fp) ) fprintf(stderr, "Failed to close %s", out_file_name);
}

void vtk_append_double(const char *out_file_name, const char *var_name, int start_index,
                       int end_index, double *values) {
    int i;
    FILE *fp = NULL;

    if ( (fp = fopen(out_file_name, "a")) == NULL ) {
        fprintf(stderr, "Failed to open %s", out_file_name);
        return;
    }

    /*
     * The first line gives the name of the dataset (variable) and its type (here "DOUBLE")
     * The second line selects the color table to use and it is usually "LOOKUP_TABLE default".
     * The following lines contain a value of the dataset per line
     */
    fprintf(fp, "SCALARS %s DOUBLE\n", var_name);
    fprintf(fp, "LOOKUP_TABLE default\n");
    for ( i = start_index; i <= end_index; i++ )
        fprintf(fp, "%f\n", values[i]);

    fprintf(fp, "\n");

    if ( fclose(fp) ) fprintf(stderr, "Failed to close %s", out_file_name);
}

void vtk_append_integer(const char *out_file_name, const char *var_name, int start_index,
                        int end_index, int *values) {
    int i;
    FILE *fp = NULL;

    fp = fopen(out_file_name, "a");
    if ( fp == NULL ) {
        fprintf(stderr, "Failed to open %s", out_file_name);
        return;
    }

    /*
     * The first line gives the name of the dataset (variable) and its type (here "INT")
     * The second line selects the color table to use and it is usually "LOOKUP_TABLE default".
     * The following lines contain a value of the dataset per line
     */
    fprintf(fp, "SCALARS %s INT\n", var_name);
    fprintf(fp, "LOOKUP_TABLE default\n");
    for ( i = start_index; i <= end_index; i++ )
        fprintf(fp, "%d\n", values[i]);

    fprintf(fp, "\n");

    if ( fclose(fp) ) fprintf(stderr, "Failed to close %s", out_file_name);
}

void visualize_distribution(char *input_file_name, char *output_name, int my_rank, int data_size, int *local_global_index)
{
    int i;

    // Setup VTK test file
    char vtk_test_file_name[50];
    char rank_num[10];
    strcpy(vtk_test_file_name, output_name);
    snprintf(rank_num, sizeof(rank_num), "_%d", my_rank);
    strcat(vtk_test_file_name, rank_num);
    strcat(vtk_test_file_name, ".vtk");

    // Setup test output to VTK
    double *test_output;
    test_output = malloc(data_size * sizeof(*test_output));
    for(i = 0; i < data_size; i++)
    {
        test_output[i] = (double)my_rank + 1.0;
    }

    test_distribution(input_file_name, vtk_test_file_name, local_global_index, data_size, test_output);

    free(test_output);
}

visualize_communication(char *input_file_name, char *output_name, char *part_type, char *comm_type, int my_rank, int nghb_idx,
    int send_cnt, int *send_lst)
{
    int i;
    double *test_output;

    // Setup VTK test file
    char vtk_test_file_name[50];
    char my_rank_num[25];
    char nghb_rank_idx[25];
    char comm[25];
    char part[25];

    strcpy(vtk_test_file_name, output_name);

    snprintf(comm, sizeof(comm), ".%s", comm_type);
    strcat(vtk_test_file_name, comm);

    snprintf(my_rank_num, sizeof(my_rank_num), ".rank%d", my_rank);
    strcat(vtk_test_file_name, my_rank_num);

    snprintf(nghb_rank_idx, sizeof(nghb_rank_idx), ".neighbor%d", nghb_idx);
    strcat(vtk_test_file_name, nghb_rank_idx);

    snprintf(part, sizeof(part), ".%s", part_type);
    strcat(vtk_test_file_name, part);

    strcat(vtk_test_file_name, ".vtk");

    // Setup test output to VTK
    test_output = calloc(send_cnt, sizeof(*test_output));
    for(i = 0; i < send_cnt; i++)
    {
        test_output[i] = (double)nghb_idx + 1.0;
    }

    test_distribution(input_file_name, vtk_test_file_name, send_lst, send_cnt, test_output);
    // printf("  %s\n", vtk_test_file_name);

    free(test_output);
}
