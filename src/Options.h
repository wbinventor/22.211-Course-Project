/*
 * Options.h
 *
 *  Created on: Jan 21, 2012
 *      Author: Lulu Li
 *				MIT, Course 22
 *              lululi@mit.edu
 *
 *  Stores global program options
 *
 */

#ifndef OPTIONS_H
#define OPTIONS_H

#include <string.h>
#include <stdlib.h>
#include "log.h"

class Options {
private:
	int _num_bins;
	int _num_batches;
	int _num_neutrons;
	int _num_threads;
	bool _test;
	bool _use_implicit_capture;
	bool _use_forced_collision;
	float _weight_low;
	float _weight_avg;
    const char* _verbosity;
public:
    Options(int argc, const char **argv);
    ~Options(void);
    int getNumBins() const;
    int getNumBatches() const;
    int getNumNeutrons() const;
    int getNumThreads() const;
    bool executeTestProblem() const;
    bool useImplicitCapture() const;
    bool useForcedCollision() const;
    float getWeightLow() const;
    float getWeightAvg() const;
    const char* getVerbosity() const;
};

#endif
