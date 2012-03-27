/*
 * Options.cpp
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  Stores global program options
 *
 */

#include "Options.h"

#define LAST(str) (strcmp(argv[i-1], (str)) == 0)


/**
 * Options constructor
 * @param argc the number of command line arguments from console
 * @param argv a char array of command line arguments from console
 */
Options::Options(int argc, const char **argv) {

	_num_bins = 100;
	_num_batches = 10;
	_num_neutrons = 1000;
	_num_threads = 1;
	_test = false;
	_use_implicit_capture = false;
	_use_forced_collision = false;
	_weight_low = 0.1;
	_weight_avg = 0.5;
	_verbosity = "NORMAL";

	for (int i = 0; i < argc; i++) {
		if (i > 0) {
			if (LAST("-nbins"))
				_num_bins = atoi(argv[i]);
			else if (LAST("-nbatches"))
				_num_batches = atoi(argv[i]);
			else if (LAST("-nneutrons"))
				_num_neutrons = atoi(argv[i]);
			else if (LAST("-nthreads"))
				_num_threads = atoi(argv[i]);
			else if (strcmp(argv[i], "-t") == 0)
				_test = true;
			else if (strcmp(argv[i], "-ic") == 0)
				_use_implicit_capture = true;
			else if (strcmp(argv[i], "-fc") == 0)
				_use_forced_collision = true;
			else if (LAST("-wavg"))
				_weight_avg = atof(argv[i]);
			else if (LAST("-wlow"))
				_weight_low = atof(argv[i]);
			else if (LAST("--verbosity") || LAST("-v"))
				_verbosity = strdup(argv[i]);
		}
	}
}

Options::~Options(void) { }


int Options::getNumBins() const {
	return _num_bins;
}


int Options::getNumBatches() const {
	return _num_batches;
}


int Options::getNumNeutrons() const {
	return _num_neutrons;
}

int Options::getNumThreads() const {
	return _num_threads;
}


bool Options::executeTestProblem() const {
	return _test;
}

bool Options::useImplicitCapture() const {
	return _use_implicit_capture;
}


bool Options::useForcedCollision() const {
	return _use_forced_collision;
}

float Options::getWeightLow() const {
	return _weight_low;
}

float Options::getWeightAvg() const {
	return _weight_avg;
}

const char* Options::getVerbosity() const {
	return _verbosity;
}
