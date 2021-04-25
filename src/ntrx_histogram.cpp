#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ntrx_histogram.h"
#include "ntrx_annot_region_tools.h"
#include "ntrx_file_utils.h"

#include <math.h>

#include <vector>
#include <algorithm>

using namespace std;

bool __DUMP_HISTOGRAM_MSGS__ = false;


void get_stats(double* vals, int n_pts, double& mean, double& var)
{
	// Get the mean.
	double total = 0.0;
	for(int i = 0; i < n_pts; i++)
	{
		total += vals[i];
	} // i loop.

	mean = total / n_pts;

	// Get the variance.
	var = 0.0;
	for(int i = 0; i < n_pts; i++)
	{
		var += (vals[i] - mean) * (vals[i] - mean);
	} // i loop.

	var = var / (n_pts - 1);
}

void get_stats(vector<double>* energies, double& mean, double& std_dev)
{
	mean = 0.0;
	for(int i = 0; i < (int)energies->size(); i++)
	{
		mean += energies->at(i);
	} // i loop.

	mean /= (int)energies->size();

	std_dev = 0.0;
	for(int i = 0; i < (int)energies->size(); i++)
	{
		std_dev += (energies->at(i) - mean) * (energies->at(i) - mean);
	} // i loop.

	// Do unbiased computation for sample variance estimate.
	std_dev /= ((int)energies->size() - 1);

	std_dev = pow(std_dev, 0.5);
}
