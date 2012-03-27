/*
 * BatchBinSet.h
 *
 *  Created on: Mar 17, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef BATCHBINSET_H_
#define BATCHBINSET_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Binner.h"
#include "log.h"


/**
 * This class represents a collection of Binner class objects. Each binner
 * is used for a specific batch of tallies for a monte carlo run. After
 * monte carlo, this class can compute batch statistics (mu, variance,
 * std dev, and rel error) and output these statistics to a file
 */
class BatchBinSet
{
private:
	char* _name;
	Binner* _binners;
	int _num_batches;
	int _num_bins;
	bool _statistics_compute;
	float* _batch_mu;
	float* _batch_variance;
	float* _batch_std_dev;
	float* _batch_rel_err;
public:
	BatchBinSet();
	virtual ~BatchBinSet();

	char* getBatchBinSetName();
	Binner* getBinner(int batch);
	float* getBatchMu();
	float* getBatchVariance();
	float* getBatchStdDev();
	float* getBatchRelativeError();

	void setBatchBinSetName(char* name);
	void createBinners(float start, float end, int num_bins,
			int num_batches, binType bin_type, tallyType tally_type,
															char* isotopes);
	void createBinners(float* bin_edges, int num_bins, int num_batches,
										tallyType tally_type, char* isotopes);

	void computeBatchStatistics();
	void computeScaledBatchStatistics(float scale_factor);
	void outputBatchStatistics(const char* filename);
};

#endif /* BATCHBINSET_H_ */
