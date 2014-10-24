/**
 * Helper functions for writing results to VTK and text files
 *
 * @author E. Xue, V. Petkov, A. Berariu
 * @date 22-May-2009, 22-Oct-2012, 22-Oct-2014
 */
#include <stdio.h>
#include <stdlib.h>
#include "util_write_files.h"



/**
 * export the results to a file
 *
 * @param infilename
 * @param outfilename
 * @param nintci
 * @param nintcf
 * @param var
 * @param iter
 * @param ratio
 * @return
 */
int write_result(char *infilename, char *outfilename, int nintci, int nintcf,
        double *var, int iter, double ratio)
{
    double *ipoint = (double *) malloc((nintcf + 1) * sizeof(double));
    int i1, i2, i3, i4, i5;

    int nc;
    for (nc = nintci; nc <= nintcf; nc++)
        ipoint[nc] = nc;

    if (nintcf <= 1)
    {
        printf("error: nintcf <= 1\n");
        return -1;
    }

    i1 = nintcf + 1;

    while (i1 != 0)
    {
        i1 = i1 / 2;
        i2 = nintcf + 1 - i1;
        i4 = 1;

        do
        {
            i3 = i4;
            do
            {
                i5 = i3 + i1;
                if (var[i3] <= var[i5])
                    break;

                double  zdum = var[i3],
                                idum  = ipoint[i3];

                var[i3] = var[i5];
                ipoint[i3] = ipoint[i5];
                var[i5] = zdum;
                ipoint[i5] = idum;
                i3 = i3 - i1;
            } while (i3 >= 1);
            i4++;
        } while (i4 < i2);
    }

    FILE *fp = fopen(outfilename, "w");
    if (fp == NULL)
    {
        printf("error opening file %s for writing\n", outfilename);
        return -1;
    }

    printf("Writing data to %s: ", outfilename);
    fprintf(fp, "               ----------------------------------------\n");
    fprintf(fp, "               - AVL -  Linear Equation Solver - GCCG -\n");
    fprintf(fp, "               ----------------------------------------\n");
    fprintf(fp, "                     diagonal scaling\n");
    fprintf(fp, "\n\n");
    fprintf(fp, "     Input File:  %s\n", infilename);
    fprintf(fp, "     ===========\n\n");
    fprintf(fp, "     Output File:  %s\n", outfilename);
    fprintf(fp, "     ============\n\n");
    fprintf(fp, "     No. of Active Cells:  %d\n", nintcf);
    fprintf(fp, "     ====================\n\n");
    fprintf(fp, "     Iteration Count - Residual Ratio\n");
    fprintf(fp, "     ================================\n\n");
    fprintf(fp, "%d %e\n", iter, ratio);
    fprintf(fp, "     Addresses Solution (Minima)           Addresses Solution (Maxima)\n");
    fprintf(fp, "     ===========================           ===========================\n\n");

    int n;
    for (n = 1; n <= 10; n++)
        fprintf(fp, "%lf %lf %lf %lf\n", ipoint[n], var[n], ipoint[nintcf - n + 1], var[nintcf - n + 1]);

    fclose(fp);
    printf("done!\n");

    free(ipoint);
    return 0;
}


int write_result_vtk(char *outFileName, int startIdx, int endIdx, int nodeCnt, int **points, int **elems, double *vector)
{
    int i,j;
    int cellCnt = endIdx - startIdx + 1;
    FILE *fp = fopen(outFileName, "w");
    if(fp == NULL)
    {
        printf("Error opening file %s for writing\n", outFileName);
        return -1;
    }

    fprintf(fp, "# vtk DataFile Version 3.0\n"); // file version and identifier
    fprintf(fp,"vtk output\n");                  // header
    fprintf(fp,"ASCII\n");                       // file format
    fprintf(fp,"DATASET UNSTRUCTURED_GRID\n");   // dataset structure

    fprintf(fp,"POINTS %d float\n",nodeCnt);

    printf("nodeCnt = %d, startInd= %d, endInd = %d, cellCnt=%d\n",nodeCnt, startIdx, endIdx, cellCnt);

    for (i = 1; i<=nodeCnt; i++){
        fprintf(fp,"%d %d %d\n",points[0][i],points[1][i],points[2][i]);
    }

    fprintf(fp,"\nCELLS %d %d\n",cellCnt, 9*cellCnt);
    for (i = startIdx; i<= endIdx; i++){
        fprintf(fp,"8 %d %d %d %d %d %d %d %d\n",elems[0][i]-1,elems[1][i]-1,
            elems[3][i]-1,elems[2][i]-1, elems[4][i]-1,elems[5][i]-1, elems[7][i]-1, elems[6][i]-1);
    }

    fprintf(fp,"\nCELL_TYPES %d\n",cellCnt);
    for (i = startIdx; i<= endIdx; i++){
        fprintf(fp,"11\n");//cell type vtk voxel
    }

    fprintf(fp,"\nCELL_DATA %d\n",cellCnt);
    fprintf(fp,"SCALARS scalars float\nLOOKUP_TABLE default\n");

    for (i = startIdx; i<= endIdx; i++){
        fprintf(fp,"%f ",vector[i]);
    }

    fclose(fp);
    return 0;
}
