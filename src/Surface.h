/*
 * Surface.h
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <vector>
#include "Region1D.h"
#include "Neutron.h"

/* Pre-define the Region1D class so the compiler knows it exists */
class Region1D;

/* The types of boundaries */
typedef enum boundaryTypes {
	REFLECTIVE,
	VACUUM,
	INTERFACE
} boundaryType;


/**
 * The Surface class represents a plane in two dimensions that is perpendicular
 * to the x-axis. The Surface knows the Region1D class objects which border
 * it and it knows whether it has reflective or vacuum boundary conditions, or
 * whether it is an interface between two regions. The Surface contains a vector
 * of neutrons which are on the Surface, and it either adds them to a bordering
 * Region1D or kills them depending on what type of boundary it is
 */
class Surface
{
private:
	float _x;
	boundaryType _boundary_type;
	Region1D* _left_region;
	Region1D* _right_region;
	std::vector<neutron*> _neutrons;
public:
	Surface();
	virtual ~Surface();

	float getX();
    boundaryType getBoundaryType() const;
    Region1D *getRightRegion() const;
    Region1D *getLeftRegion() const;

    void setX(float x);
    void setBoundaryType(boundaryType type);
    void setLeftRegion(Region1D* left_region);
    void setRightRegion(Region1D* right_region);

    void addNeutron(neutron* neutron);
    bool onSurface(float x);
    void moveNeutrons();
};

#endif /* SURFACE_H_ */
