/**
 * Helper functions for reading from input data file
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#include <stdio.h>
#include <stdlib.h>
#include "util_read_files.h"

/**
 * Parse an binary input data set and initialize simulation variables
 *
 * @param file_name
 * @param NINTCI
 * @param NINTCF
 * @param NEXTCI
 * @param NEXTCF
 * @param LCC
 * @param BS
 * @param BE
 * @param BN
 * @param BW
 * @param BL
 * @param BH
 * @param BP
 * @param SU
 * @param points_count
 * @param points
 * @param elems
 * @return
 */
int read_binary_geo(char *file_name, int *NINTCI, int *NINTCF, int *NEXTCI, int *NEXTCF, int ***LCC,
                    double **BS, double **BE, double **BN, double **BW, double **BL, double **BH,
                    double **BP, double **SU, int* points_count, int*** points, int** elems) {
    int i = 0;
    FILE *fp = fopen(file_name, "rb");

    if ( fp == NULL ) {
        fprintf(stderr, "Error opening file %s\n", file_name);
        return -1;
    }

    // 4 variables in total!!!
    fread(NINTCI, sizeof(int), 1, fp);
    fread(NINTCF, sizeof(int), 1, fp);
    fread(NEXTCI, sizeof(int), 1, fp);
    fread(NEXTCF, sizeof(int), 1, fp);

    // allocating LCC
    if ( (*LCC = (int**) malloc((*NINTCF + 1) * sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate first dimension of LCC");
        return -1;
    }

    for ( i = 0; i < *NINTCF + 1; i++ ) {
        if ( ((*LCC)[i] = (int *) malloc(6 * sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc failed to allocate second dimension of LCC\n");
            return -1;
        }
    }

    // start reading LCC
    // Note that C array index starts from 0 while Fortran starts from 1!
    for ( i = (*NINTCI); i <= *NINTCF; i++ ) {
        fread(&(*LCC)[i][0], sizeof(int), 1, fp);
        fread(&(*LCC)[i][1], sizeof(int), 1, fp);
        fread(&(*LCC)[i][2], sizeof(int), 1, fp);
        fread(&(*LCC)[i][3], sizeof(int), 1, fp);
        fread(&(*LCC)[i][4], sizeof(int), 1, fp);
        fread(&(*LCC)[i][5], sizeof(int), 1, fp);
    }

    // allocate other arrays
    if ( (*BS = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BS) failed\n");
        return -1;
    }

    if ( (*BE = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BE) failed\n");
        return -1;
    }

    if ( (*BN = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BN) failed\n");
        return -1;
    }

    if ( (*BW = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BW) failed\n");
        return -1;
    }

    if ( (*BL = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BL) failed\n");
        return -1;
    }

    if ( (*BH = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BH) failed\n");
        return -1;
    }

    if ( (*BP = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(BP) failed\n");
        return -1;
    }

    if ( (*SU = (double *) malloc((*NEXTCF + 1) * sizeof(double))) == NULL ) {
        fprintf(stderr, "malloc(SU) failed\n");
        return -1;
    }

    // read the other arrays
    for ( i = (*NINTCI); i <= *NINTCF; i++ ) {
        fread(&((*BS)[i]), sizeof(double), 1, fp);
        fread(&((*BE)[i]), sizeof(double), 1, fp);
        fread(&((*BN)[i]), sizeof(double), 1, fp);
        fread(&((*BW)[i]), sizeof(double), 1, fp);
        fread(&((*BL)[i]), sizeof(double), 1, fp);
        fread(&((*BH)[i]), sizeof(double), 1, fp);
        fread(&((*BP)[i]), sizeof(double), 1, fp);
        fread(&((*SU)[i]), sizeof(double), 1, fp);
    }

    // read geometry
    // allocate elems
    if ( (*elems = (int*) malloc((*NINTCF + 1) * 8 * sizeof(int))) == NULL ) {
        fprintf(stderr, "malloc failed to allocate elems");
        return -1;
    }

    // read elems
    for ( i = (*NINTCI); i < (*NINTCF + 1) * 8; i++ ) {
        fread(&((*elems)[i]), sizeof(int), 1, fp);
    }

    fread(points_count, sizeof(int), 1, fp);

    // allocate points vec
    if ( (*points = (int **) calloc(*points_count, sizeof(int*))) == NULL ) {
        fprintf(stderr, "malloc() POINTS 1st dim. failed\n");
        return -1;
    }

    for ( i = 0; i < *points_count; i++ ) {
        if ( ((*points)[i] = (int *) calloc(3, sizeof(int))) == NULL ) {
            fprintf(stderr, "malloc() POINTS 2nd dim. failed\n");
            return -1;
        }
    }

    int coordIdx;
    int pointIdx;
    for ( pointIdx = 0; pointIdx < *points_count; pointIdx++ ) {
        for ( coordIdx = 0; coordIdx < 3; coordIdx++ ) {
            fread(&((*points)[pointIdx][coordIdx]), sizeof(int), 1, fp);
        }
    }

    fclose(fp);

    return 0;
}


