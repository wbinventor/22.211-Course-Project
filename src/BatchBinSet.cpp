/*
 * BatchBinSet.cpp
 *
 *  Created on: Mar 17, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "BatchBinSet.h"


/**
 * Default constructor for a BatchBinSet
 */
BatchBinSet::BatchBinSet() {
	/* Default for empty batch bin set */
	_name = (char*)"";
	_num_batches = 0;
	_statistics_compute = false;
}


/**
 * BatchBinSet destructor deletes memory for arrays of binners,
 * and batch statistics if that memory has been allocated
 */
BatchBinSet::~BatchBinSet() {

	if (_num_batches != 0) {
		delete [] _binners;
		delete [] _batch_mu;
		delete [] _batch_variance;
		delete [] _batch_std_dev;
		delete [] _batch_rel_err;
	}
}


/**
 * Returns the user specified name for this BatchBinSet
 * @return the name of this BatchBinSet
 */
char* BatchBinSet::getBatchBinSetName() {
	return _name;
}


/**
 * Returns a Binner class pointer for a particular batch
 * @param batch the batch number
 * @return a pointer to a Binner class
 */
Binner* BatchBinSet::getBinner(int batch) {

	if (_num_batches == 0)
		log_printf(ERROR, "Unable to return binner %d from BatchBinSet %s "
				"since the binners have not yet been created", batch, _name);

	else if (batch < 0 || batch >= _num_batches)
		log_printf(ERROR, "Unable to return binner %d from BatchBinSet %s "
				"since this BatchBinSet has %d batches", batch, _name,
															_num_batches);

	return &_binners[batch];
}


/**
 * Returns a pointer to an array of batch averages if they have been
 * computed
 * @return a double array of batch averages for each bin
 */
double* BatchBinSet::getBatchMu() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for "
				"BatchBinSet %s so batch mu cannot be returned", _name);

	return _batch_mu;
}


/**
 * Returns a pointer to an array of batch variances if they have been
 * computed
 * @return a double array of batch variances for each bin
 */
double* BatchBinSet::getBatchVariance() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for "
				"BatchBinSet %s so batch variance cannot be returned", _name);

	return _batch_variance;
}


/**
 * Returns a pointer to an array of batch standard deviations if they have
 * been computed
 * @return a double array of batch standard deviations for each bin
 */
double* BatchBinSet::getBatchStdDev() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for "
				"BatchBinSet %s so batch std dev cannot be returned", _name);

	return _batch_std_dev;
}


/**
 * Returns a pointer to an array of batch relative errors if they have been
 * computed
 * @return a double array of batch relative errors for each bin
 */
double* BatchBinSet::getBatchRelativeError() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for "
		"BatchBinSet %s so batch relative error cannot be returned", _name);

	return _batch_rel_err;
}


/**
 * Sets this BatchBinSet's name
 * @param name a character array representing the name
 */
void BatchBinSet::setBatchBinSetName(char* name) {

	_name = name;

	/* Set the name of each Binner */
	for (int b=0; b < _num_batches; b++)
		_binners[b].setBinnerName(name);

	return;
}


/**
 * Creates a set of binners each with the same set of bins between a
 * certain start and end point
 * @param start starting bin value
 * @param end ending bin value
 * @param num_bins the number of bins for each Binner
 * @param num_batches the number of batches (number of Binner classes)
 * @param bin_type the type of bin (EQUAL, LOGARITHMIC, OTHER)
 * @param tally_type the type of tally (FLUX, ABSORB_RATE, COLLISION_RATE)
 * @param isotopes the name of the material these bins apply to or "all"
 * if they apply to all materials
 */
void BatchBinSet::createBinners(float start, float end, int num_bins,
								int num_batches, binType bin_type,
								tallyType tally_type, char* isotopes) {

	_num_batches = num_batches;
	_num_bins = num_bins;

	/* Allocate memory for binners and batch statistics arrays */
	_binners = new Binner[_num_batches];
	_batch_mu = new double[_num_bins];
	_batch_variance = new double[_num_bins];
	_batch_std_dev = new double[_num_bins];
	_batch_rel_err = new double[_num_bins];

	/* Create a Binner for each batch */
	for (int b=0; b < num_batches; b++) {
		_binners[b].generateBinEdges(start, end, _num_bins, bin_type);
		_binners[b].setTallyType(tally_type);
		_binners[b].setIsotopes(isotopes);
	}

	return;
}


/**
 * Creates a set of binners each with the same set of bins defined by
 * an array of bin edges
 * @param bin_edges array of edges which define the bins
 * @param num_bins the number of bins for each Binner
 * @param num_batches the number of batches (number of Binner classes)
 * @param tally_type the type of tally (FLUX, ABSORB_RATE, COLLISION_RATE)
 * @param isotopes the name of the material these bins apply to or "all"
 * if they apply to all materials
 */
void BatchBinSet::createBinners(float* bin_edges, int num_bins,
							int num_batches, tallyType tally_type,
												char* isotopes) {

	_num_batches = num_batches;
	_num_bins = num_bins;

	/* Allocate memory for binners and batch statistics arrays */
	_binners = new Binner[_num_batches];
	_batch_mu = new double[_num_bins];
	_batch_variance = new double[_num_bins];
	_batch_std_dev = new double[_num_bins];
	_batch_rel_err = new double[_num_bins];

	/* Create a Binner for each batch */
	for (int b=0; b < num_batches; b++) {
		_binners[b].setBinEdges(bin_edges, _num_bins);
		_binners[b].setTallyType(tally_type);
		_binners[b].setIsotopes(isotopes);
	}

	return;
}


/**
 * Computes average, variance, standard deviation and relative error for each
 * bin over the set of batches
 */
void BatchBinSet::computeBatchStatistics() {

	if (_num_batches == 0)
		log_printf(ERROR, "Cannot compute batch statistics for BatchBinSet %s"
				" since the binners have not yet been generated", _name);

	/* Loop over each bin */
	for (int i=0; i < _num_bins; i++) {

		/* Initialize statistics to zero */
		_batch_mu[i] = 0.0;
		_batch_variance[i] = 0.0;
		_batch_std_dev[i] = 0.0;
		_batch_rel_err[i] = 0.0;

		/* Accumulate flux from each batch */
		for (int j=0; j < _num_batches; j++)
			_batch_mu[i] += _binners[j].getTally(i);

		/* Compute average flux for this bin */
		_batch_mu[i] /= double(_num_batches);

		/* Compute the variance for this bin */
		for (int j=0; j < _num_batches; j++) {
			_batch_variance[i] += (_binners[j].getTally(i) - _batch_mu[i])
					* (_binners[j].getTally(i) - _batch_mu[i]);
		}
		_batch_variance[i] /= double(_num_batches);

		/* Compute the standard deviation for this bin */
		_batch_std_dev[i] = sqrt(_batch_variance[i]);

		/* Compute the relative error for this bin */
		_batch_rel_err[i] = _batch_std_dev[i] / _batch_mu[i];
	}

	_statistics_compute = true;

	return;
}


/**
 * Computes average, variance, standard deviation and relative error for each
 * bin over the set of batches. This method first scales each bin value by
 * a scaling factor
 * @param scale_factor the factor to scale each bin value by
 */
void BatchBinSet::computeScaledBatchStatistics(float scale_factor) {

	if (_num_batches == 0)
		log_printf(ERROR, "Cannot compute batch statistics for BatchBinSet %s "
				"since the binners have not yet been generated", _name);

	/* Loop over each bin */
	for (int i=0; i < _num_bins; i++) {

		/* Initialize statistics to zero */
		_batch_mu[i] = 0.0;
		_batch_variance[i] = 0.0;
		_batch_std_dev[i] = 0.0;
		_batch_rel_err[i] = 0.0;

		/* Accumulate flux from each batch */
		for (int j=0; j < _num_batches; j++)
			_batch_mu[i] += _binners[j].getTally(i) / double(scale_factor);

		/* Compute average flux for this bin */
		_batch_mu[i] /= double(_num_batches);

		/* Compute the variance for this bin */
		for (int j=0; j < _num_batches; j++) {
			_batch_variance[i] += (_binners[j].getTally(i) / scale_factor
			- _batch_mu[i]) * (_binners[j].getTally(i) / scale_factor
												- _batch_mu[i]);
		}
		_batch_variance[i] /= double(_num_batches);

		/* Compute the standard deviation for this bin */
		_batch_std_dev[i] = sqrt(_batch_variance[i]);

		/* Compute the relative error for this bin */
		_batch_rel_err[i] = _batch_std_dev[i] / _batch_mu[i];
	}

	_statistics_compute = true;

	return;
}


/**
 * Outputs the batch statistics (if they have been computed) to an
 * ASCII file
 * @param filename the output filename
 */
void BatchBinSet::outputBatchStatistics(const char* filename) {

	if (_num_batches == 0)
		log_printf(ERROR, "Cannot output batch statistics for BatchBinSet %s "
				"since the binners have not yet been generated", _name);

	else if (!_statistics_compute)
		log_printf(ERROR, "Cannot output batch statistics for BatchBinSet %s "
				"since they have not yet been computed", _name);

	/* Create output file */
	FILE* output_file;
	output_file = fopen(filename, "w");

	/* Print header to output file */
	fprintf(output_file, "Bin center, Mu, Variance, Std Dev, Rel Err\n");

	/* Loop over each bin and print mu, var, std dev and rel err */
	for (int i=0; i < _num_bins; i++) {
		fprintf(output_file, "%1.10f, %1.10f, %1.10f, %1.10f, %1.10f\n",
				_binners[0].getBinCenters()[i], _batch_mu[i],
				_batch_variance[i], _batch_std_dev[i], _batch_rel_err[i]);
	}

	fclose(output_file);

	return;
}
