/**
 * Helper functions for reading from input data file
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#include <stdio.h>
#include <stdlib.h>
#include "util_read_files.h"

int read_formatted(char *filename, char* file_format, int *nintci, int *nintcf, int *nextci,
        int *nextcf, int ***lcc, double **bs, double **be, double **bn,
        double **bw, double **bl, double **bh, double **bp, double **su)
{
    int i;
    FILE *fp;

    // Determine file format from which to read
    if(!strcmp(file_format, "text"))
    {
        fp = fopen(filename, "r");
    }
    else
    {
        fp = fopen(filename, "rb");
    }

    if (fp == NULL)
    {
        printf("Error opening file %s\n", filename);
        return -1;
    }

    //4 variables in total!!!
    if(!strcmp(file_format, "text"))
    {
        fscanf(fp, "%d", nintci);
        fscanf(fp, "%d", nintcf);
        fscanf(fp, "%d", nextci);
        fscanf(fp, "%d", nextcf);
    }
    else
    {
        fread(nintci, sizeof(int), 1, fp);
        fread(nintcf, sizeof(int), 1, fp);
        fread(nextci, sizeof(int), 1, fp);
        fread(nextcf, sizeof(int), 1, fp);
    }
    
    //allocating lcc
    if ((*lcc = (int**) malloc( ( (*nintcf) - (*nintci) + 1) * sizeof(int*) )) == NULL)
    {
        printf("malloc failed to allocate first dimension of lcc (nintcf)");
        return -1;
    }
    for (i = (*nintci); i <= (*nintcf); i++)
    {
        if (((*lcc)[i] = (int *) malloc( 6 * sizeof(int) )) == NULL)
        {
            printf("malloc(lcc) failed\n");
            return -1;
        }
    }

    //start reading lcc
    //note that c array index starts from 0 while fortran starts from 1!
    if(!strcmp(file_format, "text"))
    {
        for (i = (*nintci); i <= (*nintcf); i++)
        {
            fscanf(fp, "%d", &(*lcc)[i][0]);
            fscanf(fp, "%d", &(*lcc)[i][1]);
            fscanf(fp, "%d", &(*lcc)[i][2]);
            fscanf(fp, "%d", &(*lcc)[i][3]);
            fscanf(fp, "%d", &(*lcc)[i][4]);
            fscanf(fp, "%d", &(*lcc)[i][5]);
        }
    }
    else
    {
        for (i = (*nintci); i <= (*nintcf); i++)
        {
            fread(&(*lcc)[i][0], sizeof(int), 1, fp);
            fread(&(*lcc)[i][1], sizeof(int), 1, fp);
            fread(&(*lcc)[i][2], sizeof(int), 1, fp);
            fread(&(*lcc)[i][3], sizeof(int), 1, fp);
            fread(&(*lcc)[i][4], sizeof(int), 1, fp);
            fread(&(*lcc)[i][5], sizeof(int), 1, fp);
        }
    }

    // allocate other arrays
    if ((*bs = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*be = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bn = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bw = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bl = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bh = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*bp = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }
    if ((*su = (double *) malloc((*nextcf + 1) * sizeof(double))) == NULL)
    {
        printf("malloc() failed\n");
        return -1;
    }

    // read the other arrays
    if(!strcmp(file_format, "text"))
    {
        for (i = (*nintci); i <= *nintcf; i++)
        {
            fscanf(fp, "%lf", &((*bs)[i]));
            fscanf(fp, "%lf", &((*be)[i]));
            fscanf(fp, "%lf", &((*bn)[i]));
            fscanf(fp, "%lf", &((*bw)[i]));
            fscanf(fp, "%lf", &((*bl)[i]));
            fscanf(fp, "%lf", &((*bh)[i]));
            fscanf(fp, "%lf", &((*bp)[i]));
            fscanf(fp, "%lf", &((*su)[i]));
        }
    }
    else
    {
        for (i = (*nintci); i <= *nintcf; i++)
        {
            fread(&((*bs)[i]), sizeof(double), 1, fp);
            fread(&((*be)[i]), sizeof(double), 1, fp);
            fread(&((*bn)[i]), sizeof(double), 1, fp);
            fread(&((*bw)[i]), sizeof(double), 1, fp);
            fread(&((*bl)[i]), sizeof(double), 1, fp);
            fread(&((*bh)[i]), sizeof(double), 1, fp);
            fread(&((*bp)[i]), sizeof(double), 1, fp);
            fread(&((*su)[i]), sizeof(double), 1, fp);
        }
    }

    fclose(fp);

    return 0;
}
