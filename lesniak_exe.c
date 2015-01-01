#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdlib.h>
#include <time.h>

//#define REFRACTORY_PERIOD 1e-3 // Refractory period in s

double resistance, capacitance, voltage_threshold; // LIF parameters
double dt; // Sampling period
static const double **current_array;
size_t current_length;


// Declaring test functions
void main();
void run_code();

// Declaration of used functions
void get_spike_trace_array_lesniak(const double resistance_input, const double capacitance_input, const double voltage_threshold_input,
    const double **current_array_input, const size_t current_length_input, const double dt_input, int *spike_trace_array, 
    const int mcnc_group_no, const int *mcnc_grouping);
double rk4(const double y_old, const size_t time_index, const double (*dydt)(double, double, int), const int mcnc_group_id);
double get_d_voltage_d_time(const double time_index, const double voltage, const int mcnc_group_id);
double gaussian_noise(const double mean, const double std);


void main()
{
	clock_t begin, end;
	double time_spent;
	begin = clock();
	run_code();
	end = clock();
	time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
	printf("Total time: %f sec.", time_spent);
	return;
}


void run_code()
{
	FILE *file_ptr;
	double **current_array_input;
	int *spike_trace_array;
	double resistance_input, capacitance_input, voltage_threshold_input;
	double buffer;
	double dt_input = 1e-4;
    int mcnc_group_no = 4;
    int mcnc_grouping[4] = {8, 5, 3, 1};
	size_t i, current_length_input;
	// Initialize model parameters
	resistance_input = 5e8;
	capacitance_input = 1e-11;
	voltage_threshold_input = 29e-3;
	// Initialize arrays
	current_length_input = 50001;
	spike_trace_array = (int *)calloc(current_length_input, sizeof(int));
    current_array_input = (double **)malloc(sizeof(double *) * current_length_input);
    for(i=0; i<current_length_input; i++)
        current_array_input[i] = (double *)calloc(sizeof(double), mcnc_group_no);
	// Read the file into current_input
	file_ptr = fopen("../test/csvs/trans_current.csv", "r");
	for(i=0; i<current_length_input; i++)
	{
		fscanf(file_ptr, "%lf,%lf,%lf,%lf,%lf\n", &(buffer), &(current_array_input[i][0]), &(current_array_input[i][1]),
            &(current_array_input[i][2]), &(current_array_input[i][3]));
	}
	fclose(file_ptr);
	get_spike_trace_array_lesniak(resistance_input, capacitance_input, voltage_threshold_input, current_array_input, current_length_input, dt_input, 
        spike_trace_array, mcnc_group_no, mcnc_grouping);
	for(i=0; i<current_length_input; i++)
		if(spike_trace_array[i] == 1)
			printf("%f\n", i*dt_input);
	free(spike_trace_array);
    for(i=0; i<current_length_input; i++)
        free(current_array_input[i]);
    free(current_array_input);
}


void get_spike_trace_array_lesniak(const double resistance_input, const double capacitance_input, const double voltage_threshold_input,
    const double **current_array_input, const size_t current_length_input, const double dt_input, int *spike_trace_array, 
    const int mcnc_group_no, const int *mcnc_grouping)
{
	/* 
	The wrapper for the lif model. 
	1. The refractory period is abstracted out - because in real data, almost no ISIs are below 1 ms. 
	2. The output needs to be initialized to all zeros before passing in!
	*/
	size_t i;
    int mcnc_group_id, fire;
	double *voltage_new, *voltage_old;
	const double (*dydt)(double, double, int);
	// Assign the value for model parameters (global variables)
	resistance = resistance_input;
	capacitance = capacitance_input;
	voltage_threshold = voltage_threshold_input;
	// Assign value to sampling period
	dt = dt_input;
	// Assign the value for the current_array input
	current_array = current_array_input;
	current_length = current_length_input;
	// Assign the governing equation
	dydt = &get_d_voltage_d_time;
	// Get the output value
	voltage_old = (double *)calloc(mcnc_group_no, sizeof(double));
    voltage_new = (double *)calloc(mcnc_group_no, sizeof(double));
	fire = 0;
	for(i=1; i<current_length-1; i++)
	{
        // Calculate new voltage for all branches and check firing
        for(mcnc_group_id=0; mcnc_group_id<mcnc_group_no; mcnc_group_id++)
        {
            voltage_new[mcnc_group_id] = rk4(voltage_old[mcnc_group_id], i, dydt, mcnc_group_id);
            if(voltage_new[mcnc_group_id] > voltage_threshold)
            {
                fire = 1;
                break;
            }
        }
        // If fired, then reset all branches and record timing
		if(fire)
		{
			spike_trace_array[i] = 1;
			fire = 0;
            for(mcnc_group_id=0; mcnc_group_id<mcnc_group_no; mcnc_group_id++)
            {
                voltage_old[mcnc_group_id] = 0.;
            }
		}
		else
		{
            for(mcnc_group_id=0; mcnc_group_id<mcnc_group_no; mcnc_group_id++)
            {
                voltage_old[mcnc_group_id] = voltage_new[mcnc_group_id];
            }
		}
	}
	return;
}


double rk4(const double y_old, const size_t time_index, const double (*dydt)(double, double, int), const int mcnc_group_id)
{
	/* 
	ODE stepper using Runge-Kutta 4th order method. The input boundary condition needs to be uniformly sampled.

	Parameters
	----------
	double y_old:
		The y value at the beginning of the step.
	size_t time_index:
		The index of the y_old value, used to know what time it is at the beginning of the step.
	double (*dydt)(double, double):
		dydt is a function pointer, pointing to dy/dt. The 1st input parameter is the time index, which could be either integer or integer + 0.5; 
		2nd input parameter is simply y.

	Return
	------
	double y_new:
		The y value at the end of the step.
	*/
	double k1, k2, k3, k4; // RK4 intermediate parameters
	double y_new;
	k1 = (*dydt)(time_index+0., y_old, mcnc_group_id);
	k2 = (*dydt)(time_index+.5, y_old+dt/2.*k1, mcnc_group_id);
	k3 = (*dydt)(time_index+.5, y_old+dt/2.*k2, mcnc_group_id);
	k4 = (*dydt)(time_index+1., y_old+dt*k3, mcnc_group_id);
	y_new = y_old + 1./6. * dt * (k1 + 2 * k2 + 2 * k3 + k4);
	return y_new;
}


double get_d_voltage_d_time(const double time_index, const double voltage, const int mcnc_group_id)
{
	double current, d_voltage_d_time;
	current = .5 * current_array[(size_t)floor(time_index)][mcnc_group_id] + .5 * current_array[(size_t)ceil(time_index)][mcnc_group_id];
	d_voltage_d_time = -1. / (resistance * capacitance) * voltage + current / capacitance;
	return d_voltage_d_time;
}


double gaussian_noise(const double mean, const double std)
{
    static int have_spare = 0;
    static double u1, u2, z1, z2;
    if(have_spare)
    {
        have_spare = 0;
        z2 = sqrt(-2. * log(u1)) * sin(2. * M_PI * u2);
        return mean + std * z2;
    }
    have_spare = 1;
    u1 = ((double) (rand() + 1) / (RAND_MAX + 1));
    u2 = ((double) (rand() + 1) / (RAND_MAX + 1));
    z1 = sqrt(-2. * log(u1)) * cos(2. * M_PI * u2);
    return mean + std * z1;
}