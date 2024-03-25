#include <math.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <assert.h>
#include <string.h>
#include "tools.h"

#define SQ(X) ((X) * (X))
#define CB(X) ((X) * (X) * (X))
#define QD(X) ((X) * (X) * (X) * (X))
#define PI 3.1415926536
#define N_DIMS 3


double wave_function(double *r1, double *r2, double alpha)
{
	double r1_norm = vector_norm(r1, N_DIMS);
	double r2_norm = vector_norm(r2, N_DIMS);

	double r12 = distance_between_vectors(r1, r2, N_DIMS);

	return exp(-2. * (r1_norm + r2_norm) + r12 / ( 2.0 * (1.0 + alpha  * r12)));
}


double local_energy(double *r1, double *r2, double alpha)
{
	double r1_unit[N_DIMS], r2_unit[N_DIMS];

	memcpy(r1_unit, r1, sizeof(double[N_DIMS]));
	memcpy(r2_unit, r2, sizeof(double[N_DIMS]));

	normalize_vector(r1_unit, N_DIMS);
	normalize_vector(r2_unit, N_DIMS);

	double r12 = distance_between_vectors(r1, r2, N_DIMS);

	double dot_product = 0.0;

	for(int j = 0; j < N_DIMS; ++j){
		dot_product += (r1_unit[j] - r2_unit[j]) * (r1[j] - r2[j]);
	}

	return -4.0 + dot_product / (r12 * SQ(1.0 + alpha * r12))
					   	- 1.0 / (r12 * CB(1.0 + alpha * r12))
					   	- 1.0 / (4.  * QD(1.0 + alpha * r12))
					   	+ 1.0 / r12;
}


void local_energy_series(double *energy, double positions[][2 * N_DIMS], double alpha, int N)
{
	double r1[N_DIMS], r2[N_DIMS];
	int i, j;

	for (i = 0; i < N; ++i){
		for(j = 0; j < N_DIMS; ++j){
			r1[j] = positions[i][j];
			r2[j] = positions[i][j+N_DIMS];
		}
		energy[i] = local_energy(r1, r2, alpha);
	}
}


void metropolis_method(double positions[][2 * N_DIMS], double d, double alpha, int N)
{
	const gsl_rng_type *T;
	const int seed = 42;

	gsl_rng *r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, seed);

	const double init_displacement = 50.0;

	double r1_new[N_DIMS], 	r1[N_DIMS];
	double r2_new[N_DIMS],	r2[N_DIMS];
	double p1, p2;
	double accept_ratio = 0.;
	double d_step, xi;
	int accepts = 0;

	int j, k, m, n;
	for (j = 0; j < N_DIMS; ++j){
		r1[j] 		= init_displacement * (gsl_rng_uniform(r) - 0.5);
		r2[j] 		= init_displacement * (gsl_rng_uniform(r) - 0.5);
		r1_new[j] 	= r1[j];
		r2_new[j] 	= r2[j];
	}

	for (int i = 0; i < N; ++i){
		d_step = d * (gsl_rng_uniform(r) - 0.5);

		m = (int) (gsl_rng_uniform(r) * 2);
		n = (int) (gsl_rng_uniform(r) * N_DIMS);

		if (m == 1)
			r1_new[n] = r1[n] + d_step;
		else
			r2_new[n] = r2[n] + d_step;

		p1 = SQ(wave_function(r1, 		r2, 	alpha));
		p2 = SQ(wave_function(r1_new, 	r2_new, alpha));

		accept_ratio = p2 / p1;
		xi = gsl_rng_uniform(r);
		accepts += (accept_ratio >= xi);

		if (accept_ratio >= xi){
			for (k = 0; k < N_DIMS; ++k){
				r1[k] = r1_new[k];
				r2[k] = r2_new[k];
			}
		}
		else {
			for (k = 0; k < N_DIMS; ++k){
				r1_new[k] = r1[k];
				r2_new[k] = r2[k];
			}
		}

		for (k = 0; k < N_DIMS; ++k){
			positions[i][k] 	= r1[k];
			positions[i][k + N_DIMS] = r2[k];
		}
	}

	//printf("Acceptance ratio %f\n", (double) accepts / N);

	gsl_rng_free(r);
}


void correlation_theta(double *corr_array, double positions[][2 * N_DIMS], int N)
{
	double r1[N_DIMS], 	r2[N_DIMS];
	double r1_norm, r2_norm;
	int i, j;

	for (i = 0; i < N; ++i){
		for (j = 0; j < N_DIMS; ++j){
			r1[j] = positions[i][j];
			r2[j] = positions[i][j + N_DIMS];
		}

		r1_norm = vector_norm(r1, N_DIMS);
		r2_norm = vector_norm(r2, N_DIMS);

		corr_array[i] = dot_product(r1, r2, N_DIMS) / (r1_norm * r2_norm);
	}
}


double correlation_function(double *arr, double *corr_arr, int M_c, int N)
{
	double corr_mean;
	int i, j;

	double correlation 	= 0.;
	double sum 			= 0.;

	double var 			= SQ(standard_deviation(arr, N));
	double mean 		= average(arr, N);

	for (i = 0; i < M_c; ++i){
		corr_mean = 0;
		for (j = 0; j < N; ++j){
			corr_mean += arr[j] * arr[j + i];
		}
		corr_mean /= N;
		correlation = (corr_mean - SQ(mean))/ var;
		sum += correlation;
		corr_arr[i] = correlation;
	}
	return 2.0 * sum;
}


void block_averaging(double *n_s, double *energy, int N, int B_max)
{
	int i, j, M_B;
	double var_averages, var_energies;
	double *averages = calloc(N, sizeof(double));


	for (int B = 1; B < B_max; ++B){
		M_B = N / B;
		for (i = 0; i < M_B; ++i){
			averages[i] = 0.;
			for (j = 0; j < B; ++j){
				averages[i] += energy[j + i * B] / B;
			}
		}
		var_averages 	= SQ(standard_deviation(averages, M_B));
		var_energies 	= SQ(standard_deviation(energy, N));
		n_s[B-1] 		= (B * var_averages) / var_energies;
	}
	free(averages); 	averages = NULL;
}


void remove_eq_data(double *array, int N_eq, int N)
{
	double *temp_array = malloc(sizeof(double[N]));

	memcpy(temp_array, array, sizeof(double[N]));

	for (int i = 0; i < N - N_eq; ++i){
		array[i] = temp_array[i + N_eq];
	}
	memset(array + (N - N_eq), 0, sizeof(double[N_eq]));

	free(temp_array);	temp_array = NULL;
}


void alpha_sweep(int n_cycles, int N_eq, int n_alpha, double n_s, double d)
{

	int n_rel = n_cycles - N_eq;
	double alpha0 = 0.05, alpha1 = 0.4;
	double alpha[n_alpha];
	double mean_array[n_alpha];
	double std_array[n_alpha];

	linspace(alpha, alpha0, alpha1, n_alpha);

	double (*positions)[2 * N_DIMS]	= malloc(sizeof(double[n_cycles][2 * N_DIMS]));
	double *local_energies 			= malloc(sizeof(double[n_cycles]));

	for (int i = 0; i < n_alpha; ++i){
		memset(positions, 0, sizeof(double[n_cycles][2 * N_DIMS]));
		memset(local_energies, 0, n_cycles * sizeof(double));

		metropolis_method(positions, d, alpha[i], n_cycles);

		local_energy_series(local_energies, positions, alpha[i], n_cycles);

		remove_eq_data(local_energies, N_eq, n_cycles);

		mean_array[i] = average(local_energies, n_rel);
		std_array[i]  = standard_deviation(local_energies, n_rel) / sqrt(n_rel / n_s);
		printf("%d\n",i);
	}

	double **f_write_array = malloc(sizeof(double[N_DIMS][n_rel]));

	f_write_array[0] = alpha;
	f_write_array[1] = mean_array;
	f_write_array[2] = std_array;

	write_to_file("data/alpha_sweep.txt",f_write_array, n_alpha, N_DIMS);

	free(f_write_array);	f_write_array	= NULL;
	free(local_energies); 	local_energies 	= NULL;
	free(positions);	  	positions		= NULL;
}


void alpha_grad_descent(int n_cycles, int N_eq, int N_max, double beta, double d)
{
	int n_rel = n_cycles - N_eq;
	double alpha[N_max + 1];
	alpha[0] = 0.16;
	alpha[1] = 0.15;

	double energy_local;
	double log_gradient, log_wavefn1, log_wavefn2;
	double sum_gradient, sum_energy,  sum_energy_wavefn;

	int j,k;
	double r1[N_DIMS], r2[N_DIMS];

	double (*pos)[2 * N_DIMS] = malloc(sizeof(double[n_cycles][2 * N_DIMS]));

	for (int i = 0; i < N_max; ++i){
		memset(pos, 0, sizeof(double[n_cycles][2 * N_DIMS]));
		metropolis_method(pos, d, alpha[i], n_cycles);
		if (i >= 1){
			sum_gradient		= 0.;
			sum_energy 			= 0.;
			sum_energy_wavefn	= 0.;

			for (j = N_eq; j < n_cycles; ++j){

				for (k = 0; k < N_DIMS; ++k){
					r1[k] = pos[j][k];
					r2[k] = pos[j][k + N_DIMS];
				}
				energy_local 		 = local_energy(r1, r2, alpha[i]);

				log_wavefn1  		 = log(wave_function(r1, r2, alpha[i-1]));
				log_wavefn2  		 = log(wave_function(r1, r2, alpha[i]));

				log_gradient 		 = (log_wavefn2 -  log_wavefn1) / (alpha[i] - alpha[i-1]);
				sum_gradient 		+= log_gradient;

				sum_energy 			+= energy_local;
				sum_energy_wavefn 	+= energy_local * log_gradient;
			}
			sum_energy_wavefn 	/= n_rel;
			sum_energy			/= n_rel;
			sum_gradient 		/= n_rel;

			double grad_alpha = 2.0 * (sum_energy_wavefn - sum_energy * sum_gradient);
			alpha[i + 1] = alpha[i] - pow(i, -beta) * grad_alpha;

			assert(alpha[i + 1] > 0.0);
		}
	}

	double **f_alpha = malloc((N_max + 1) * sizeof(double));
	f_alpha[0] = alpha;

	write_to_file("data/grad_descent_0_10.txt", f_alpha, N_max + 1, 1);

	free(f_alpha);	f_alpha = NULL;
	free(pos); 		pos 	= NULL;
}





