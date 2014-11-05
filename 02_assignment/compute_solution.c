/**
 * Computational loop
 *
 * @file compute_solution.c
 * @date 22-Oct-2012
 * @author V. Petkov
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int compute_solution(const int max_iters, int nintci, int nintcf, int nextcf, int** lcc, double* bp,
                     double* bs, double* bw, double* bl, double* bn, double* be, double* bh,
                     double* cnorm, double* var, double *su, double* cgup, double* residual_ratio,
                     int* local_global_index, int* global_local_index, int neighbors_count,
                     int* send_count, int** send_list, int* recv_count, int** recv_list) {
    /** parameters used in gccg */
    int iter = 1;
    int if1 = 0;
    int if2 = 0;
    int nor = 1;
    int nor1 = nor - 1;
    int nc = 0;
    int nomax = 3;

    /** the reference residual */
    double resref = 0.0;

    /** array storing residuals */
    double *resvec = (double *) calloc(sizeof(double), (nintcf + 1));

    // initialize the reference residual
    for ( nc = nintci; nc <= nintcf; nc++ ) {
        resvec[nc] = su[nc];
        resref = resref + resvec[nc] * resvec[nc];
    }

    resref = sqrt(resref);
    if ( resref < 1.0e-15 ) {
        fprintf(stderr, "Residue sum less than 1.e-15 - %lf\n", resref);
        return 0;
    }

    /** the computation vectors */
    double *direc1 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *direc2 = (double *) calloc(sizeof(double), (nextcf + 1));
    double *adxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *adxor2 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor1 = (double *) calloc(sizeof(double), (nintcf + 1));
    double *dxor2 = (double *) calloc(sizeof(double), (nintcf + 1));

    while ( iter < max_iters ) {
        /**********  START COMP PHASE 1 **********/
        // update the old values of direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc1[nc] = direc1[nc] + resvec[nc] * cgup[nc];
        }

        // compute new guess (approximation) for direc
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            direc2[nc] = bp[nc] * direc1[nc] - bs[nc] * direc1[lcc[nc][0]]
                         - be[nc] * direc1[lcc[nc][1]] - bn[nc] * direc1[lcc[nc][2]]
                         - bw[nc] * direc1[lcc[nc][3]] - bl[nc] * direc1[lcc[nc][4]]
                         - bh[nc] * direc1[lcc[nc][5]];
        }
        /********** END COMP PHASE 1 **********/

        /********** START COMP PHASE 2 **********/
        // execute normalization steps
        double oc1, oc2, occ;
        if ( nor1 == 1 ) {
            oc1 = 0;
            occ = 0;

            for ( nc = nintci; nc <= nintcf; nc++ ) {
                occ = occ + direc2[nc] * adxor1[nc];
            }

            oc1 = occ / cnorm[1];
            for ( nc = nintci; nc <= nintcf; nc++ ) {
                direc2[nc] = direc2[nc] - oc1 * adxor1[nc];
                direc1[nc] = direc1[nc] - oc1 * dxor1[nc];
            }

            if1++;
        } else {
            if ( nor1 == 2 ) {
                oc1 = 0;
                occ = 0;

                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + direc2[nc] * adxor1[nc];
                }

                oc1 = occ / cnorm[1];
                oc2 = 0;
                occ = 0;
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    occ = occ + direc2[nc] * adxor2[nc];
                }

                oc2 = occ / cnorm[2];
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    direc1[nc] = direc1[nc] - oc1 * dxor1[nc] - oc2 * dxor2[nc];
                    direc2[nc] = direc2[nc] - oc1 * adxor1[nc] - oc2 * adxor2[nc];
                }

                if2++;
            }
        }

        // compute the new residual
        cnorm[nor] = 0;
        double omega = 0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            cnorm[nor] = cnorm[nor] + direc2[nc] * direc2[nc];
            omega = omega + resvec[nc] * direc2[nc];
        }

        omega = omega / cnorm[nor];
        double res_updated = 0.0;
        for ( nc = nintci; nc <= nintcf; nc++ ) {
            resvec[nc] = resvec[nc] - omega * direc2[nc];
            res_updated = res_updated + resvec[nc] * resvec[nc];
            var[nc] = var[nc] + omega * direc1[nc];
        }

        res_updated = sqrt(res_updated);
        *residual_ratio = res_updated / resref;

        // exit on no improvements of residual
        if ( *residual_ratio <= 1.0e-10 ) break;

        iter++;

        // prepare additional arrays for the next iteration step
        if ( nor == nomax ) {
            nor = 1;
        } else {
            if ( nor == 1 ) {
                for ( nc = nintci; nc <= nintcf; nc++ ) {
                    dxor1[nc] = direc1[nc];
                    adxor1[nc] = direc2[nc];
                }
            } else {
                if ( nor == 2 ) {
                    for ( nc = nintci; nc <= nintcf; nc++ ) {
                        dxor2[nc] = direc1[nc];
                        adxor2[nc] = direc2[nc];
                    }
                }
            }

            nor++;
        }
        nor1 = nor - 1;
        /********** END COMP PHASE 2 **********/
    }

    free(direc1);
    free(direc2);
    free(adxor1);
    free(adxor2);
    free(dxor1);
    free(dxor2);
    free(resvec);

    return iter;
}


