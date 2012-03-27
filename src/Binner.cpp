/*
 * binner.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Binner.h"


/**
 * Default Binner constructor
 */
Binner::Binner() {

	_name = (char*)"";

	 /* Sets the default delta between bins to zero */
	_bin_delta = 0;

	/* Default is to tally all isotopes */
	_isotopes = (char*)"all";
}


/**
 * Binner destructor deletes memory for tallies, number of tallies,
 * bin centers and bin edges if they have been created
 */
Binner::~Binner() {

	if (_num_bins != 0) {
		delete [] _tallies;
		delete [] _num_tallies;
		delete [] _centers;
		if (_bin_type != OTHER)
			delete [] _edges;
	}
}


/**
 * Returns the name of this Binner as specified by the user
 * @return the Binner's name
 */
char* Binner::getBinnerName() {
	return _name;
}


/**
 * Returns the number of bins
 * @return the number of bins
 */
int Binner::getNumBins() {
	return _num_bins;
}


/**
 * Returns a float array of bin edge values
 * @return array of bin edge values
 */
float* Binner::getBinEdges() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return bin edges for Binner %s since "
				 "the bins have not yet been created", _name);

	 return _edges;
}


/**
 * Returns a float array of bin center values
 * @return array of bin center values
 */
float* Binner::getBinCenters() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return bin centers for Binner %s since the "
				 "bins have not yet been created", _name);

	 return _centers;
}


/**
 * Returns the delta spacing between bins. NOTE: this value is only non-zero
 * for EQUAL and LOGARITHMIC bin types
 * @return the spacing between bins
 */
float Binner::getBinDelta() {
	return _bin_delta;
}


float Binner::getBinDelta(float sample) {
	/* If this Binner uses equally spaced bins in linear or logarithmic
	 * space, return the bin delta */
	if (_bin_type == EQUAL || _bin_type == LOGARITHMIC)
		return _bin_delta;

	/* If instead this Binner uses irregularly spaced bin edges defined
	 * by a users, compute bin delta of the bin around the sample */
	else {
		int bin_index = getBinIndex(sample);
		return (_edges[bin_index] - _edges[bin_index-1]);
	}
}


/**
 * Returns the type of bin (EQUAL, LOGARITHMIC, OTHER)
 * @return the bin type
 */
binType Binner::getBinType() {
	return _bin_type;
}


/**
 * Returns the type of tally for these bins (FLUX, COLLISION_RATE,
 * or ABSORB_RATE)
 * @return the tally type
 */
tallyType Binner::getTallyType() {
	return _tally_type;
}


/**
 * Returns a float array of the tallies within each bin
 * @return an array of
 */
float* Binner::getTallies() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return tallies for Binner %s since the "
				 "bins for have not yet been created", _name);

	 return _tallies;
}


/**
 * Returns a specific tally for a specific bin
 * @param bin_index the index for the bin of interest
 * @return the tally within that bin
 */
float Binner::getTally(int bin_index) {

	if (bin_index < 0 || bin_index >= _num_bins)
		log_printf(ERROR, "Tried to get a tally for a bin index for Binner %s"
				"which does not exist: %d, num_bins = %d", _name, bin_index,
				_num_bins);

	return _tallies[bin_index];
}


/**
 * Returns an int array of the number of times tallied within each bin
 * @return an array of the number of tallies in each bin
 */
int* Binner::getNumTallies() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return tally numbers for Binner %s since "
				 "the bins have not yet been created", _name);

	return _num_tallies;
}


/**
 * Returns the number of times tallied within a specific bin
 * @param bin_index the bin of interest
 * @return the number of tallies in that bin
 */
int Binner::getNumTallies(int bin_index) {

	if (bin_index < 0 || bin_index >= _num_bins)
		log_printf(ERROR, "Tried to get a tally number for Binner %s for "
				"a bin index which does not exist: %d", _name, bin_index);

	return _num_tallies[bin_index];
}


/**
 * Returns the maximum tally value over all bins
 * @return the maximum tally value
 */
float Binner::getMaxTally() {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return the maximum tally for Binner %s"
				 "since the bins have not yet been created", _name);

	float max_tally = 0;

	/* Loop over all bins */
	for (int i=0; i < _num_bins; i++) {
		if (_tallies[i] > max_tally)
			max_tally = _tallies[i];
	}

	return max_tally;
}


/**
 * Returns the maximum tally value over all bins
 * @return the maximum tally value
 */
float Binner::getMinTally() {
	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return the minimum tally for Binner %s"
				 " since the bins have not yet been created", _name);

	float min_tally = std::numeric_limits<float>::infinity();

	/* Loop over all bins */
	for (int i=0; i < _num_bins; i++) {
		if (_tallies[i] < min_tally)
			min_tally = _tallies[i];
	}

	return min_tally;
}


/**
 * Finds the bin index for a sample in a set of bins. If the samples
 * is outside the bounds of all bins, it returns infinity
 * @param sample the sample value of interest
 * @return the bin index for the sample
 */
int Binner::getBinIndex(float sample) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return a bin index for Binner %s since "
				 "the bins have not yet been created", _name);

	/* Set index to infinity to begin with */
	int index = std::numeric_limits<float>::infinity();

	/* if the sample is equal to the last bin edge, return the last bin */
	if (sample == _edges[_num_bins])
		return _num_bins-1;

	/* Equally spaced bins */
	if (_bin_type == EQUAL)
		index = int((sample - _edges[0]) / _bin_delta);

	/* Logarithmically spaced bins */
	else if (_bin_type == LOGARITHMIC)
		index = int((log10(sample) - log10(_edges[0])) / _bin_delta);

	/* If the bin_type == OTHER then the bin edges were not generated by
	 * generateEqualBinEdges, so use a brute force search to find the bin */
	else {

		/* Loop over all bin edges to find the correct bin index */
		for (int i=0; i <= _num_bins; i++) {
			if (sample >= _edges[i] && sample < _edges[i+1]) {
				index = i;
				break;
			}
		}
	}

	/* If this sample was not contained within a bin set index to infinity*/
	if (index > _num_bins)
		index = std::numeric_limits<float>::infinity();

	return index;
}


/**
 * Return the isotopes that this binner is meant to tally
 * @return a character array of the isotope's name or "all" for all isotopes
 */
char* Binner::getIsotopes() {
	return _isotopes;
}



/**
 * Sets this Binner's name
 * @param name the name of the Binner
 */
void Binner::setBinnerName(char* name) {
	_name = name;
}


/**
 * Set the type of tally for this Binner (FLUX, COLLISION_RATE, ABSORB_RATE)
 * @param type the tally type
 */
void Binner::setTallyType(tallyType type) {
	_tally_type = type;
}


/**
 * Set a user-defined float array of bin edge values
 * @param edges the array of bin edges
 * @param num_bins the number of bins
 */
void Binner::setBinEdges(float* edges, int num_bins) {

	_num_bins = num_bins;
	_edges = edges;
	_bin_type = OTHER;

	/* Set all tallies to zero by default */
	_tallies = new float[_num_bins];
	_num_tallies = new int[_num_bins];

	/* Loop over tallies and set to zero */
	for (int i=0; i < _num_bins; i++) {
		_tallies[i] = 0;
		_num_tallies[i] = 0;
	}

	/* Create an array of the center values between bins */
	generateBinCenters();
}



/**
 * Set the isotope that this binner is meant to tally
 * @param isotopes a character array of the isotope's name or
 * "all" for all isotopes
 */
void Binner::setIsotopes(char* isotopes) {
	_isotopes = isotopes;
}


/**
 * Generate edges between bins defined by a start and end point
 * @param start first bin edge value
 * @param end last bin edge value
 * @param num_bins the number of bins to be created
 * @param type the type of bins (EQUAL or LOGARITHMIC)
 */
void Binner::generateBinEdges(float start, float end, int num_bins,
														binType type) {
	if (start == end)
		log_printf(ERROR, "Unable to create bins for Binner %s between"
				"the same start and end points: %f", _name, start);

	_num_bins = num_bins;
	_bin_type = type;

	/* Allocate memory for tallies */
	_tallies = new float[num_bins];
	_num_tallies = new int[num_bins];

	/* Set all tallies to zero by default */
	for (int i=0; i < num_bins; i++) {
		_tallies[i] = 0;
		_num_tallies[i] = 0;
	}

	/* Equal spacing between bins */
	if (type == EQUAL) {
		_bin_delta = float(end - start) / float(_num_bins);

		/* Generate points from start to end for each bin edge */
		_edges = linspace(start, end, num_bins+1);
	}

	/* Logarithmically equal spacing between bins */
	else if (type == LOGARITHMIC) {
		_bin_delta = float(log10(end) - log10(start)) / float(_num_bins);

		/* Generate points from start to end for each bin edge */
		_edges = logspace(start, end, num_bins+1);
	}

	else
		log_printf(ERROR, "Bin type %d is not yet implemented for Binner %s",
															_name, type);

	/* Create an array of the center values between bins */
	generateBinCenters();

	return;
}


/**
 * Compute the center points between bin edges for this Binner's bins
 */
void Binner::generateBinCenters() {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot generate bin centers for Binner %s since "
				 "the bins have not yet been created", _name);

	/* Allocate memory for the bin centers array */
	_centers = new float[_num_bins];

	/* Loop over all bins and find the midpoint between edges */
	for (int i=0; i < _num_bins; i++)
		_centers[i] = (_edges[i] + _edges[i+1]) / 2.0;

	return;
}


/**
 * Tallies unity for each sample in a float array of samples
 * @param samples array of samples to tally
 * @param num_samples the number of samples to tally
 */
void Binner::tally(float* samples, int num_samples) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally samples in Binner %s since the "
				 "bins have not yet been created", _name);

	int bin_index;

	/* Loop over and tally all samples */
	for (int i=0; i < num_samples; i++) {
		bin_index = getBinIndex(samples[i]);
		if (bin_index >= 0 && bin_index < _num_bins) {
			_tallies[bin_index]++;
			_num_tallies[bin_index]++;
		}
	}

	return;
}


/**
 * Tallies unity for a sample
 * @param samples array of samples to tally
 */
void Binner::tally(float sample) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(sample);

	if (bin_index >= 0 && bin_index < _num_bins) {
		_tallies[bin_index]++;
		_num_tallies[bin_index]++;
	}

	return;
}


/**
 * Tallies a weight for each sample in a float array of samples
 * @param samples array of samples to tally
 * @param sample_weights array of sample weights to increment tallies by
 * @param num_samples the number of samples to tally
 */
void Binner::weightedTally(float* samples, float* sample_weights,
														int num_samples) {
	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted samples in Binner %s "
				 "since the bins have not yet been created", _name);

	int bin_index;

	/* Loop over and tally all samples */
	for (int i=0; i < num_samples; i++) {
		bin_index = getBinIndex(samples[i]);
		if (bin_index >= 0 && bin_index < _num_bins) {
			_tallies[bin_index] += sample_weights[i];
			_num_tallies[bin_index]++;
		}
	}

	return;
}


/**
 * Tallies a weight for a sample
 * @param sample a sample to tally
 * @param weight the weight to increment tally by
 */
void Binner::weightedTally(float sample, float weight) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(sample);

	if (bin_index >= 0 && bin_index < _num_bins) {
		_tallies[bin_index] += weight;
		_num_tallies[bin_index]++;
	}

	return;
}


/**
 * Divide each tally by the maximum tally value
 */
void Binner::normalizeTallies() {

	if (_num_bins == 0)
		log_printf(ERROR, "Cannot normalize tallies for Binner %s since it is"
						"the bins have not yet been created", _name);

	float max_tally = getMaxTally();

	/* Divide each tally by maximum tally value */
	for (int n=0; n < _num_bins; n++)
		_tallies[n] /= max_tally;

	return;
}


/**
 * Divide each tally by a given scaling factor
 * @param scale_factor factor to normalize tallies by
 */
void Binner::normalizeTallies(float scale_factor) {

	if (_num_bins == 0)
		log_printf(ERROR, "Cannot normalize tallies for Binner %s since it is"
						"the bins have not yet been created", _name);

	/* Divide each tally by maximum tally value */
	for (int n=0; n < _num_bins; n++)
		_tallies[n] /= scale_factor;

	return;
}


/**
 * Helper function to generate an array of equally spaced floats between
 * a given start and end point. Modeled after MATLAB's linspace function
 * @param start the starting point
 * @param end the ending point
 * @param num_values the number of values to create
 * @return a pointer to the array of points
 */
float* linspace(float start, float end, int num_values) {

	float* values = new float[num_values];

	/* Spacing between values */
	float delta = float(end - start) / float(num_values-1);

	/* Loop over all values */
	for (int i=0; i <= num_values; i++)
		values[i] = delta * i + start;

	return values;
}


/**
 * Helper function to generate an array of equal logarithmically spaced
 * floats between a given start and end point. Modeled after MATLAB's
 * linspace function
 * @param start the starting point
 * @param end the ending point
 * @param num_values the number of values to create
 * @return a pointer to the array of points
 */
float* logspace(float start, float end, int num_values) {

	/* Create an equally spaced array of base 10 exponent values */
	float* values = linspace((float)log10(start), log10(end), num_values);

	/* Loop over all values and project back original space */
	for (int i=0; i < num_values; i++)
		values[i] = pow(10, values[i]);

	return values;
}
