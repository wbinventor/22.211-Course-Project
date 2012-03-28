/*
 * binner.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef BINNER_H_
#define BINNER_H_

#include <limits>
#include <math.h>
#include "log.h"
#include "arraycreator.h"


/* Bin spacing types */
typedef enum binTypes {
	EQUAL,
	LOGARITHMIC,
	OTHER
} binType;


/* Type of tallies */
typedef enum tallyTypes {
	FLUX_SPATIAL,
	FLUX_ENERGY,
	CAPTURE_RATE_SPATIAL,
	CAPTURE_RATE_ENERGY,
	ABSORPTION_RATE_SPATIAL,
	ABSORPTION_RATE_ENERGY,
	ELASTIC_RATE_SPATIAL,
	ELASTIC_RATE_ENERGY,
	FISSION_RATE_SPATIAL,
	FISSION_RATE_ENERGY,
	TRANSPORT_RATE_SPATIAL,
	TRANSPORT_RATE_ENERGY,
	COLLISION_RATE_SPATIAL,
	COLLISION_RATE_ENERGY,
	DIFFUSION_RATE_SPATIAL,
	DIFFUSION_RATE_ENERGY
} tallyType;


/**
 * This class represents a set of bins which are defined by a set of values
 * defining the edges between bins. This class holds the edges, the centers
 * between bins. It also allows for tallies to be made within each bin.
 */
class Binner{
private:
	char* _name;
	int _num_bins;
	float* _edges;
	double* _centers;
	double* _tallies;
	int* _num_tallies;
	float _bin_delta;
	binType _bin_type;
	tallyType _tally_type;
	char* _isotopes;
public:
	Binner();
	virtual ~Binner();
	char* getBinnerName();
	int getNumBins();
	float* getBinEdges();
	double* getBinCenters();
	float getBinDelta();
	float getBinDelta(float sample);
	binType getBinType();
	tallyType getTallyType();
	double* getTallies();
	double getTally(int bin_index);
	int* getNumTallies();
	int getNumTallies(int bin_index);
	double getMaxTally();
	double getMinTally();
	int getBinIndex(float sample);
	char* getIsotopes();

	void setBinnerName(char* name);
	void setTallyType(tallyType type);
	void setBinEdges(float* edges, int num_edges);
	void setIsotopes(char* isotopes);

	void generateBinEdges(float start, float end, int num_bins, binType type);
	void generateBinCenters();
	void tally(float* samples, int num_samples);
	void tally(float sample);
	void weightedTally(float* samples, float* sample_weights, int num_samples);
	void weightedTally(float sample, float weight);
	void normalizeTallies();
	void normalizeTallies(float scale_factor);
};

#endif /* BINNER_H_ */
