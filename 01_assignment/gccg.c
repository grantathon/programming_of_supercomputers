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
	// setup PAPI varaibles
	long_long values[NUM_EVENTS];
	long_long start_cycles, end_cycles, start_usec, end_usec;

	// { FloatingPointOps, L2_TotalAccess, L3_TotalAccess, L2_TotalMisses, L3_TotalMisses }
	//unsigned int events[NUM_EVENTS] = { PAPI_FP_OPS, PAPI_L2_TCA, PAPI_L3_TCA, PAPI_L2_TCM, PAPI_L3_TCM };
	unsigned int events[NUM_EVENTS] = { PAPI_L2_TCA, PAPI_L3_TCA, PAPI_L2_TCM, PAPI_L3_TCM, PAPI_L1_DCM };

	int i;

	const int max_iters = 10000;    /// maximum number of iteration to perform

	/** Simulation parameters parsed from the input datasets */
	int nintci, nintcf;    /// internal cells start and end index
	/// external cells start and end index. The external cells are only ghost cells.
	/// They are accessed only through internal cells
	int nextci, nextcf;
	int **lcc;    /// link cell-to-cell array - stores neighboring information
	/// Boundary coefficients for each volume cell (South, East, North, West, High, Low)
	double *bs, *be, *bn, *bw, *bh, *bl;
	double *bp;    /// Pole coefficient
	double *su;    /// Source values

	double residual_ratio;    /// the ratio between the reference and the current residual
	double *var;    /// the variation vector -> keeps the result in the end

	/** Additional vectors required for the computation */
	double *cgup, *oc, *cnorm;

	  
 
	char *file_in = argv[1];
 
	/********** START INITIALIZATION **********/
	// read-in the input file
	int init_status = initialization(file_in, &nintci, &nintcf, &nextci, &nextcf, &lcc,
									 &bs, &be, &bn, &bw, &bl, &bh, &bp, &su, &var, &cgup, &oc, 
									 &cnorm);

	if ( init_status != 0 ) {
		fprintf(stderr, "\n Failed to initialize data!\n");
		abort();
	} 
   
	/********** END INITIALIZATION **********/

	//	TODO: initialize PAPI library
	if (PAPI_library_init ( PAPI_VER_CURRENT ) != PAPI_VER_CURRENT )
	{
		printf("\nerror: PAPI initialization failed! Exiting.\n");
		exit(1);
	}

	start_cycles = PAPI_get_real_cyc();	// Gets the starting time in clock cycles
	start_usec = PAPI_get_real_usec();	// Gets the starting time in microseconds

	// start the PAPI counters. Monitoring starts after this command!
	int initErr = PAPI_start_counters( (int*)events, NUM_EVENTS);
	if (initErr != PAPI_OK)
	{
		handle_error(initErr);
		//printf("\nerror: PAPI did not start the counters!\n");
	}

	/********** START COMPUTATIONAL LOOP **********/
	int total_iters = compute_solution(max_iters, nintci, nintcf, nextcf, lcc, bp, bs, bw, bl, bn,
									   be, bh, cnorm, var, su, cgup, &residual_ratio);
	/********** END COMPUTATIONAL LOOP **********/

	// stop counters and gather the collected data
	int retErr = PAPI_stop_counters(values, NUM_EVENTS);
	if (retErr != PAPI_OK)
	{
		handle_error(retErr);
		//printf("\nerror: PAPI error while stopping the counters!\n");
	}

	// get execution time
	end_cycles = PAPI_get_real_cyc();	//Gets the ending time in clock cycles
	end_usec = PAPI_get_real_usec();	//Gets the ending time in microseconds
	
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

 
	// calculate and print the PAPI performance variables
	// { L2_TotalAccess, L3_TotalAccess, L2_TotalMisses, L3_TotalMisses, FloatingPointOps }
	float missRateL2 = ( (float)values[2] / (float)values[0] );
	float missRateL3 = ( (float)values[3] / (float)values[1] );

	printf( "Wall clock time in usecs: %lld\n", end_usec - start_usec);
	printf( "Time in clock cycles: %lld\n", end_cycles - start_cycles);	
	//printf("Mflops: %f \n", ( (float)values[0] / 1000000));
	printf("L1 misses: %lli \n", values[4]);
	printf("L2 misses: %lli \n", values[2]);
	printf("L3 misses: %lli \n", values[3]);
	printf("L2 total access: %lli \n", values[0]);
	printf("L3 total access: %lli \n", values[1]);
	printf("L2 miss rate [percent]: %f \n", (missRateL2 * 100));
	printf("L3 miss rate [percent]: %f \n", (missRateL3 * 100));

	return 0;
}

