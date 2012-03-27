/*
 * Surface.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Surface.h"


/**
 * Surface constructor sets default values
 */
Surface::Surface() {
	_boundary_type = VACUUM;		/* Default boundary type */
	_left_region = NULL;
	_right_region = NULL;
}


/**
 * Surface destructor
 */
Surface::~Surface() { }


/**
 * Returns the x value at which this Surface exists
 * @param the x value for this Surface
 */
float Surface::getX() {
	return _x;
}


/**
 * Returns this Surface's boundary type (REFLECTIVE, VACUUM, or INTERFACE)
 * @return Surface boundary type
 */
boundaryType Surface::getBoundaryType() const {
    return _boundary_type;
}


/**
 * Returns a pointer to the Region1D which borders this Surface on the left
 * @return a pointer to the Region1D
 */
Region1D* Surface::getLeftRegion() const {
    return _left_region;
}


/**
 * Returns a pointer to the Region1D which borders this Surface on the right
 * @return a pointer to the Region1D
 */
Region1D* Surface::getRightRegion() const {
    return _right_region;
}


/**
 * Sets the x value for where this Surface is located
 * @param x the x value
 */
void Surface::setX(float x) {
	_x = x;
}


/**
 * Sets the boundary type for this Surface (REFLECTIVE, VACUUM, or INTERFACE)
 * @param type the boundary type
 */
void Surface::setBoundaryType(boundaryType type) {
    _boundary_type = type;
}


/**
 * Sets the Region1D which borders this Surface to the left
 * @param left_region pointer to a Region1D class object
 */
void Surface::setLeftRegion(Region1D* left_region) {
    _left_region = left_region;
}


/**
 * Sets the Region1D which borders this Surface to the right
 * @param right_region pointer to a Region1D class object
 */
void Surface::setRightRegion(Region1D* right_region) {
    _right_region = right_region;
}


/**
 * Adds a neutron to this surface
 * @param neutron a pointer to a neutron struct
 */
void Surface::addNeutron(neutron* neutron) {
	neutron->_x = _x;
	_neutrons.push_back(neutron);
}


/**
 * Checks whether a certain x values is on the Surface
 * (within numerical error)
 * @param x the value to check
 * @return if on the Surface (true), otherwise false
 */
bool Surface::onSurface(float x) {
	if (fabs(_x - x) < 1E-6)
		return true;

	return false;
}


/**
 * Moves neutrons to either the left or right bordering Region1D
 * based on the neutron's trajectory. If this Surface has VACUUM
 * boundary conditions, it instead kills the neutron
 */
void Surface::moveNeutrons() {

	log_printf(DEBUG, "Inside Surface moveNeutrons method");

	neutron* curr;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {
		curr = (*iter);

		log_printf(DEBUG, "Moving neutron with x = %f, mu = %f", curr->_x, curr->_mu);

		/* If the neutron is not on the Surface, print error message */
		if (!onSurface(curr->_x))
			log_printf(ERROR, "Cannot move a neutron with x = %f off of "
					"Surface x = %f since it is not on the Surface",
					curr->_x, _x);

		else {
			/* If the surface is a vacuum, kill neutron */
			if (_boundary_type == VACUUM) {
				iter = _neutrons.erase(iter);
				--iter;
				delete curr;
			}

			/* If the surface is reflective */
			else if (_boundary_type == REFLECTIVE) {

				/* Reverse mu for reflection */
				curr->_mu *= -1.0;

				/* Figure out which region to put neutron in */
				if (curr->_mu <= 0)
					_left_region->addNeutron(curr);
				else
					_right_region->addNeutron(curr);
			}

			/* If the surface is an interface between two regions */
			else {

				/* Figure out which region to put neutron in */
				if (curr->_mu <= 0)
					_left_region->addNeutron(curr);
				else
					_right_region->addNeutron(curr);
			}
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}
