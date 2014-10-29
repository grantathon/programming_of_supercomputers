/**
 * Main GCCG program
 *
 * @author E. Xue, V. Petkov
 * @date 22-May-2009, 22-Oct-2012
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "initialization.h"
#include "compute_solution.h"
#include "finalization.h"
#include "vol2mesh.h"
#include "util_write_files.h"

// PAPI header file
 #include "papi.h"

// number of events to profile
 #define NUM_EVENTS 5

// PAPI error handle function
void handle_error(int retval)
{
    printf("PAPI error %d: %s\n", retval, PAPI_strerror(retval));
    exit(1);
}

int main(int argc, char *argv[])
{
    /* setup PAPI varaibles     */
    long_long eventValues[NUM_EVENTS];
    long_long start_cycles, end_cycles, start_usec[2], end_usec[2];
    int EventSet = PAPI_NULL;

    /* Initialize PAPI library for any PAPI low-level calls */
    if (PAPI_library_init ( PAPI_VER_CURRENT ) != PAPI_VER_CURRENT )
    {
        fprintf(stderr, "PAPI initialization failed! Exiting\n");
        exit(1);
    }

    /* initialize a multiplex EventSet for MFLOPS */
    if (PAPI_multiplex_init() != PAPI_OK)
    {
        fprintf(stderr, "PAPI multiplex not initialized! Exiting\n");
        exit(1);
    }

    /* Create an EventSet */
    if (PAPI_create_eventset(&EventSet) != PAPI_OK )
    {
        fprintf(stderr, "PAPI event set not created! Exiting\n");        
        exit(1);
    }

    /* attach events to PAPI_assign_eventset_component */
    if (PAPI_assign_eventset_component(EventSet, 0) != PAPI_OK)
    {
        fprintf(stderr, "PAPI failed to assign event! Exiting\n");
        exit(1);
    }

    /* turn events into  multiplex Events */
    if (PAPI_set_multiplex(EventSet) != PAPI_OK)
    {
        fprintf(stderr, "PAPI failed to set event! Exiting\n");
        exit(1);
    }

    /* add the events to the PAPI event set */
    if (PAPI_add_event(EventSet, PAPI_L2_TCA ) != PAPI_OK)  exit(1);
    if (PAPI_add_event(EventSet, PAPI_L3_TCA ) != PAPI_OK)  exit(1);
    if (PAPI_add_event(EventSet, PAPI_L2_TCM ) != PAPI_OK)  exit(1);
    if (PAPI_add_event(EventSet, PAPI_L3_TCM ) != PAPI_OK)  exit(1);
    if (PAPI_add_event(EventSet, PAPI_FP_OPS ) != PAPI_OK)  exit(1);

    int i;
    const int max_iters = 10000;    /// maximum number of iteration to perform

    /** Simulation parameters parsed from the input datasets */
    int nintci, nintcf;     /* internal cells start and end index
                             * external cells start and end index. The external cells are only ghost cells.
                             * They are accessed only through internal cells
                             */
    
    int nextci, nextcf;
    int **lcc;          /// link cell-to-cell array - stores neighboring information
                        /// Boundary coefficients for each volume cell (South, East, North, West, High, Low)

    double *bs, *be, *bn, *bw, *bh, *bl;
    double *bp;    /// Pole coefficient
    double *su;    /// Source values

    double residual_ratio;      /// the ratio between the reference and the current residual
    double *var;                /// the variation vector -> keeps the result in the end

    /** Additional vectors required for the computation */
    double *cgup, *oc, *cnorm;
    int node_count;
    int **points, **elems;

    // Check if number of arguements is correct
    if(argc != 4)
    {
        fprintf(stderr, "Please provide three arguements: <format> <input file> <output prefix>.\n");
        abort();
    }

    // Read input arguements
    char *file_format = argv[1];
    char *file_in = argv[2];
    char *output_prefix = argv[3];

    // Check if file format is valid
    if(strcmp(file_format, "text") && strcmp(file_format, "bin"))
    {
        fprintf(stderr, "Please provide a file format as either 'text' or 'bin' (first arguement).\n");
        abort();
    }
 
    /* start timer for data read   */
    start_usec[0] = PAPI_get_real_usec();  

    /********** START INITIALIZATION **********/
    // read-in the input file
    int init_status = initialization(file_in, file_format, &nintci, &nintcf, &nextci, &nextcf, &lcc,
                                     &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &var, &cgup, &oc,
                                     &cnorm);

    if ( init_status != 0 ) {
        fprintf(stderr, "Failed to initialize data!\n");
        abort();
    } 
    
    /* stop timer for data read   */
    end_usec[0] = PAPI_get_real_usec();
    
    // Display to user some input file parameters
    printf("nintci = %d\n", nintci);
    printf("nintcf = %d\n", nintcf);
    printf("nextci = %d\n", nextci);
    printf("nextcf = %d\n", nextcf);

    /********** END INITIALIZATION **********/


    /* start time and cycle measurements just before Computation Phase  */
    start_cycles = PAPI_get_real_cyc(); // Gets the starting time in clock cycles
    start_usec[1] = PAPI_get_real_usec();  // Gets the starting time in microseconds

    /*  PAPI: start the counters. Monitoring starts after this command! */
    if (PAPI_start(EventSet) != PAPI_OK)
    {
        fprintf(stderr, "PAPI failed to start the counters! Exiting\n");
        exit(1);
    }

    /********** START COMPUTATIONAL LOOP **********/
    int total_iters = compute_solution(max_iters, nintci, nintcf, nextcf, lcc, bp, bs, bw, bl, bn,
                                       be, bh, cnorm, var, su, cgup, &residual_ratio);
    /********** END COMPUTATIONAL LOOP **********/

    /* get execution Time   */
    end_cycles = PAPI_get_real_cyc();   //Gets the ending time in clock cycles
    end_usec[1] = PAPI_get_real_usec();    //Gets the ending time in microseconds

    /* Read the collected data   */
    if (PAPI_read(EventSet, eventValues) != PAPI_OK)
    {
        fprintf(stderr, "PAPI failed to read the counters! Exiting\n");
        exit(1);
    }

    /* Read again and stop counting     */
    if (PAPI_stop(EventSet, eventValues) != PAPI_OK)
    {
        fprintf(stderr, "PAPI failed to sto the counters! Exiting\n");        
        exit(1);
    }

    /* Convert the volume data to mesh data for vtk visualization output */
    vol2mesh(nintci, nintcf, lcc, &node_count, &points, &elems);

    // Output vtk files with appropriate naming convention
    char var_vtk_name[50];
    char su_vtk_name[50];
    char cgup_vtk_name[50];
    sprintf(var_vtk_name, "%s.VAR.vtk", output_prefix);
    sprintf(su_vtk_name, "%s.SU.vtk", output_prefix);
    sprintf(cgup_vtk_name, "%s.CGUP.vtk", output_prefix);
    write_result_vtk(var_vtk_name, nintci, nintcf, node_count, points, elems, var);
    write_result_vtk(su_vtk_name, nintci, nintcf, node_count, points, elems, su);
    write_result_vtk(cgup_vtk_name, nintci, nintcf, node_count, points, elems, cgup);
    
    /********** START FINALIZATION **********/
    finalization(file_in, total_iters, residual_ratio, nintci, nintcf, var, cgup, su);
    /********** END FINALIZATION **********/

    free(cnorm);
    free(var);
    free(cgup);
    free(su);
    free(bp);
    free(bh);
    free(bl);
    free(bw);
    free(bn);
    free(be);
    free(bs);

    for ( i = 0; i < 6; ++i)
    {
       printf("lcc[1][%d] = %d\n", i, lcc[1][i] );
    }

    for ( i = nintci; i <= nintcf; i++ ) {
        free(lcc[i]);
    }
    free(lcc);

    /* calculate and print the PAPI performance variables   */
    float missRate_L2 = ( (float)eventValues[2] / (float)eventValues[0] );
    float missRate_L3 = ( (float)eventValues[3] / (float)eventValues[1] );

    printf("\nPerformance Data\nMflops: %lld\n", (eventValues[4] / 1000000) );
    printf("L2 misses: %lli \n", eventValues[2]);
    printf("L3 misses: %lli \n", eventValues[3]);
    printf("L2 total access: %lli \n", eventValues[0]);
    printf("L3 total access: %lli \n", eventValues[1]);
    printf("L2 miss rate [percent]: %f \n", (missRate_L2 * 100));
    printf("L3 miss rate [percent]: %f \n", (missRate_L3 * 100));
    printf("Computation time [ms]: %f\n", (float)( end_usec[1] - start_usec[1] ) / 1000);
    printf("Time to read data file [ms]: %f\n", (float)( end_usec[0] - start_usec[0] ) / 1000);
    printf("Time in clock cycles: %lld\n", end_cycles - start_cycles);  

    return 0;
}
