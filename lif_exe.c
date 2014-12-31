#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

//#define REFRACTORY_PERIOD 1e-3 // Refractory period in s

double resistance, capacitance, voltage_threshold; // LIF parameters
double dt; // Sampling period
const double *current_array; 
size_t current_length;


// Declaring test functions
void main();
void run_code();

// Declaration of used functions
void get_spike_trace_array_lif(double resistance_input, double capacitance_input, double voltage_threshold_input, const double *current_array_input, size_t current_length_input, 
	double dt_input, int *spike_trace_array);
double rk4(double y_old, size_t time_index, double (*dydt)(double, double));
double get_d_voltage_d_time(double time_index, double voltage);

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
	double *current_array_input;
	int *spike_trace_array;
	double resistance_input, capacitance_input, voltage_threshold_input;
	double buffer;
	double dt_input = 1e-4;
	size_t i, current_length_input;
	// Initialize model parameters
	resistance_input = 5e8;
	capacitance_input = 1e-11;
	voltage_threshold_input = 29e-3;
	// Initialize arrays
	current_length_input = 50001;
	current_array_input = (double *)calloc(current_length_input, sizeof(double));
	spike_trace_array = (int *)calloc(current_length_input, sizeof(int));
	// Read the file into current_input
	file_ptr = fopen("../test/csvs/test_current.csv", "r");
	for(i=0; i<current_length_input; i++)
	{
		fscanf(file_ptr, "%lf,%lf\r\n", &(buffer), &(current_array_input[i]));
		current_array_input[i] *= 1e-3;
	}
	fclose(file_ptr);
	get_spike_trace_array_lif(resistance_input, capacitance_input, voltage_threshold_input, current_array_input, current_length_input, dt_input, spike_trace_array);
	for(i=0; i<current_length_input; i++)
		if(spike_trace_array[i] == 1)
			printf("%f\n", i*dt_input);
	free(current_array_input);
	free(spike_trace_array);
}


void get_spike_trace_array_lif(double resistance_input, double capacitance_input, double voltage_threshold_input, const double *current_array_input, size_t current_length_input, 
	double dt_input, int *spike_trace_array)
{
	/* 
	The wrapper for the lif model. 
	1. The refractory period is abstracted out - because in real data, almost no ISIs are below 1 ms. 
	2. The output needs to be initialized to all zeros before passing in!
	*/
	size_t i;
	double voltage_new, voltage_old;
	double (*dydt)(double, double);
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
	voltage_old = 0.;
	for(i=1; i<current_length; i++)
	{
		voltage_new = rk4(voltage_old, i, dydt);
		if(voltage_new>voltage_threshold)
		{
			spike_trace_array[i] = 1;
			voltage_old = 0;
		}
		else
		{
			voltage_old = voltage_new;
		}
	}
	return;
}


double rk4(double y_old, size_t time_index, double (*dydt)(double, double))
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
	k1 = (*dydt)(time_index+0., y_old);
	k2 = (*dydt)(time_index+.5, y_old+dt/2.*k1);
	k3 = (*dydt)(time_index+.5, y_old+dt/2.*k2);
	k4 = (*dydt)(time_index+1., y_old+dt*k3);
	y_new = y_old + 1./6. * dt * (k1 + 2 * k2 + 2 * k3 + k4);
	return y_new;
}


double get_d_voltage_d_time(double time_index, double voltage)
{
	double current, d_voltage_d_time;
	current = .5 * current_array[(size_t)floor(time_index)] + .5 * current_array[(size_t)ceil(time_index)];
	d_voltage_d_time = -1. / (resistance * capacitance) * voltage + current / capacitance;
	return d_voltage_d_time;
}

