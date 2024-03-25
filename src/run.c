#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "metropolis.h"
#include "tools.h"


/*
 * Global variable for equilibration time and step size.
 */
const int N_eq = 1000;
const double d = 4.0;

/*
 * Main run program.
 * Runs metropolis algorithm and determines
 * statistical inefficiency with block averaging and
 * correlation function.
 *
 * Saves local energies, correlations cos \theta,
 * positions to three txt files.
 */
void run_write_data_metropolis(int block_average_bool, int corr_fn_bool)
{
	const int n_cycles 	= 15e6;
	const double alpha 	= 0.1;

   /*
	* positions:	 	matrix (N steps x (3 dims x 2 particles))
	* local_energies:	(N steps) array to store local energies
	* correlations:		(N steps) array to store the correlations of the two electrons
	* f_array: 			temporary array to write energies
	* f_positions:		temporary array to write positions
	*/

	double (*positions)[6] 	= malloc(sizeof(double[n_cycles][6]));
	double *local_energies 	= calloc(n_cycles, sizeof(double));
	double *correlations 	= calloc(n_cycles, sizeof(double));

	double **f_energies 	= malloc(sizeof(double*));
	double **f_correlations = malloc(sizeof(double*));

	/*
	 * Runs metropolis algorithm
	 */
	metropolis_method(positions, d, alpha, n_cycles);

	/*
	 * Calculates local energies for each step
	 */
	local_energy_series(local_energies, positions, alpha, n_cycles);

	/*
	 * Calculates cos \theta for each step for histogram
	 */
	correlation_theta(correlations, positions, n_cycles);

	double **f_positions = transpose_array(n_cycles, 6, positions);

	f_energies[0] 		= local_energies;
	f_correlations[0] 	= correlations;


	write_to_file("data/energies.txt", 	f_energies, 	n_cycles, 1);
	write_to_file("data/corr.txt", 		f_correlations, n_cycles, 1);
	write_to_file("data/positions.txt", f_positions, 	n_cycles, 6);

	free(f_energies);		f_energies 		= NULL;
	free(f_correlations);	f_correlations	= NULL;
	free(f_positions);		f_positions 	= NULL;

	free(positions);	  	positions		= NULL;
	free(correlations);	  	correlations	= NULL;
	/*
	 * Removes the first samples to account for equilibration time,
	 * adds the same amount of zeros to the end of the array as it removed.
	 */
	remove_eq_data(local_energies, N_eq, n_cycles);

	/*
	 * Correlation function to determine statistical inefficiency
	 */
	if (corr_fn_bool)
	{
		const int M_c 		= 80;
		const int N_samples = 3e6;

		double correlation_arr[M_c];
		double n_s = correlation_function(local_energies, correlation_arr, M_c, N_samples);
		printf("Correlation function: n_s = %f\n", n_s);

		double **f_corr_fn = malloc(M_c * sizeof(double));
		f_corr_fn[0] = correlation_arr;

		write_to_file("data/correlation_fn.txt", f_corr_fn, M_c, 1);
		free(f_corr_fn);		f_corr_fn 		= NULL;
	}


	/*
	 * Block averaging to determine statistical inefficiency
	 */
	if(block_average_bool){
		const int B_max = 2000;
		double block_array[B_max-1];
		const int n_rel = n_cycles - N_eq;

		block_averaging(block_array, local_energies, n_rel, B_max);
		printf("Block average: n_s = %f\n", block_array[B_max-2]);

		double **f_block_array = malloc((B_max-1) * sizeof(double));
		f_block_array[0] = block_array;

		write_to_file("data/blocks.txt", f_block_array, (B_max-1), 1);
		free(f_block_array); 	f_block_array 	= NULL;
	}

	/*
	 * Free memory allocations
	 */
	free(local_energies); 	local_energies 	= NULL;
}

void run_alpha_sweep()
{
	/*
	 * n_cycles: 	number of max steps in metropolis algorithm.
	 * n_alpha: 	number of points in linear space from 0.05 to 0.4 for alpha sweep
	 * n_s: 		statistical inefficiency
	 */
	const int n_cycles 	= (int) 15e6;
	const int n_alpha 	= 50;
	const double n_s 	= 16;

	alpha_sweep(n_cycles, N_eq, n_alpha, n_s, d);
}

void run_alpha_grad_descent()
{
	const int n_cycles 	= (int) 15e6;
	const int N_max		= 20;
	const double beta 	= 0.10;
	alpha_grad_descent(n_cycles, N_eq, N_max, beta, d);
}

int run(int argc, char *argv[])
{
	//run_alpha_sweep();
	run_alpha_grad_descent();

	//const int run_block 			= 1;
	//const int run_correlation_fn 	= 1;
	//run_write_data_metropolis(run_block, run_correlation_fn);

	return 0;
}




