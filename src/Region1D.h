/*
 * Region1D.h
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef REGION1D_H_
#define REGION1D_H_

#include <vector>
#include <math.h>
#include <algorithm>
#include <stdarg.h>
#include <omp.h>
#include "Material.h"
#include "Isotope.h"
#include "Neutron.h"
#include "Surface.h"
#include "Binner.h"

/* Pre-define Surface class so compiler knows it exists */
class Surface;


typedef enum varianceReduceTypes {
	IMPLICIT_CAPTURE,
	FORCED_COLLISION
} varianceReduceType;


typedef enum regionTypes {
	BOUNDED,
	INFINITE,
	TWO_REGION_PIN_CELL
} regionType;


	/**
 * The Region1D class represents a single dimensional region
 * bounded by two planes, or Surface class objects. The region
 * contains a vector of neutrons which live within it and is filled
 * by a set of Isotope class objects, each with a different number
 * density. The Region1D class contains all of the physics for moving
 * colliding neutrons
 */
class Region1D
{
private:
	regionType _region_type;
	char* _region_name;

	/* If region is BOUNDED type then it has a left and right Surface */
	Surface* _left_surface;
	Surface* _right_surface;

	/* Vector of pointers to neutrons contained by this Region1D */
	std::vector<neutron*> _neutrons;

	/* Vector of pointers to neutrons transferred to this Region1D
	 * during a two-region pin cell simulation */
	std::vector<neutron*> _transfer_neutrons;

	/* Binners for tallying */
	std::vector<Binner*> _bin_sets;

	Material* _material;
	float _volume;

	/* Two region pin cell parameters */
	bool _fuel;
	bool _moderator;
	Region1D* _other_region;
	float _sigma_e;
	float _beta;
	float _alpha1;
	float _alpha2;

	/* Variance reduction */
	bool _use_implicit_capture;
	bool _use_forced_collision;
	float _weight_low;
	float _weight_avg;
public:
	Region1D();
	virtual ~Region1D();

    Surface* getLeftSurface() const;
    Surface* getRightSurface() const;
    float getVolume();
    Material* getMaterial();
    regionType getRegionType();
    char* getRegionName();
    bool isFuel();
    bool isModerator();
    Region1D* getOtherPinCellRegion();

    void setRegionName(char* region_name);
    void setMaterial(Material* material);
    void setLeftSurface(Surface* surface);
    void setRightSurface(Surface* surface);
    void setVolume(float volume);
    void addBinner(Binner* bins);
    void setAsFuel();
    void setAsModerator();
    void setOtherPinCellRegion(Region1D* _other_region);
    void setTwoRegionPinCellParams(float sigma_e, float beta, float alpha1,
    															float alpha2);
    void useImplicitCapture(float weight_low, float weight_avg);
    void useForcedCollision(float weight_low, float weight_avg);

    bool playRussianRoulette(neutron* neutron);
    float computeFuelFuelCollisionProb(int energy_index);
    float computeModeratorFuelCollisionProb(int energy_index);
    void clearBinners();
    bool contains(float x);
    bool onBoundary(float x);
    void addNeutron(neutron* neutron);
    void transferNeutron(neutron* neutron);
    void initializeTransferredNeutrons();
    void twoRegionNeutronTransferral();
    void moveNeutrons();
    int getNumNeutrons();
};

#endif /* REGION1D_H_ */
