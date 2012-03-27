/*
 * Integrator.cpp
 *
 *  Created on: Feb 9, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "integrate.h"


/**
 * This method performs a 2D numerical integral over arrays of x and y values
 * using a particular integration method
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 * @param sceheme integration method
 * @return the result of the numerical integral
 */
float integrate(float* x, float* y, int length, integrationScheme scheme) {

	float integral = 0;

	switch(scheme) {
		case RIEMANN_RIGHT:
			integral = computeRiemannRight(x, y, length);
			break;
		case RIEMANN_LEFT:
			integral = computeRiemannLeft(x, y, length);
			break;
		case RIEMANN_CENTER:
			integral = computeRiemannCenter(x, y, length);
			break;
		case TRAPEZOIDAL:
			integral = computeTrapezoidal(x, y, length);
			break;
		case SIMPSONS:
			integral = computeSimpsons(x, y, length);
			break;
		case SIMPSONS38:
			integral = computeSimpsons38(x, y, length);
			break;
		case BOOLES:
			integral = computeBooles(x, y, length);
			break;
		case ENUM_END:
			break;
	}

	return integral;
}


/**
 * This method performs a cumulative numerical integral over arrays of x and
 * y values using a particular integration method
 * @param x the the x values
 * @param y the y values
 * @param cdf the array of cdf values at each value of x and y
 * @param length the length of the x and y arrays
 * @param sceheme integration method
 */
void cumulativeIntegral(float* x, float* y, float* cdf, int length,
													integrationScheme scheme) {

	/* Calculate cumulative integral */
	for (int i=1; i < length+1; i++)
		cdf[i-1] = integrate(x, y, i, scheme);
}



/**
 * This method performs a numerical integral using the left-centered Riemann
 * numerical integration method
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 */
float computeRiemannLeft(float* x, float* y, int length) {

	float integral = 0;
	float delta_x = 0;

	for (int i = 0; i < length; i++) {
		if (i < length - 1)
			delta_x = x[i+1] - x[i];
		else
			delta_x = x[i] - x[i-1];

		integral += delta_x * y[i];
	}

	return integral;
}


/**
 * This method performs a numerical integral using the right-centered Riemann
 * numerical integration method
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 */
float computeRiemannRight(float* x, float* y, int length) {

	float integral = 0;
	float delta_x = 0;

	for (int i = 1; i < length-1; i++) {
		delta_x = x[i] - x[i-1];

		integral += delta_x * y[i];
	}

	return integral;
}


/**
 * This method performs a numerical integral using the midpoint-centered
 * Riemann numerical integration method
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 */
float computeRiemannCenter(float* x, float* y, int length) {

	float integral = 0;
	float delta_x = 0;

	for (int i = 0; i < length; i++) {
		if (i == 0)
			delta_x = x[i+1] - x[i];
		else if (i == length - 1)
			delta_x = x[i] - x[i-1];
		else
			delta_x = (x[i] - x[i-1]) / 2.0 + (x[i+1] - x[i]) / 2.0;

		integral += delta_x * y[i];
	}

	return integral;
}


/**
 * This method performs a numerical integral using the trapezoidal numerical
 * integration method
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 */
float computeTrapezoidal(float* x, float* y, int length) {

	float integral = 0;
	float delta_x = 0;

	for (int i = 1; i < length; i++) {
		delta_x = x[i] - x[i-1];
		integral += delta_x * (y[i] + y[i-1]) / 2.0;
	}

	return integral;
}


/**
 * This method performs a numerical integral using the Simpson's numerical
 * integration method
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 */
float computeSimpsons(float* x, float* y, int length) {

	float integral = 0;
	float delta_x = 0;

	for (int i = 0; i < length - 2; i++) {
		delta_x = x[i+1] - x[i];
		integral += (delta_x / 6.0) * (y[i] + 4*y[i+1] + y[i+2]);
	}

	return integral;
}


/**
 * This method performs a numerical integral using the Simpson's 3/8ths
 * numerical integration method
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 */
float computeSimpsons38(float* x, float* y, int length) {

	float integral = 0;
	float delta_x = 0;

	for (int i = 0; i < length - 3; i++) {
		delta_x = x[i+1] - x[i];
		integral += (delta_x / 8.0) * (y[i] + 3*y[i+1] + 3*y[i+2] + y[i+3]);
	}

	return integral;
}


/**
 * This method performs a numerical integral using the Boole's numerical
 * integration method
 * @param x the the x values
 * @param y the y values
 * @param length the length of the x and y arrays
 */
float computeBooles(float* x, float* y, int length) {

	float integral = 0;
	float delta_x = 0;

	for (int i = 0; i < length - 4; i++) {
		delta_x = x[i+1] - x[i];
		integral += (delta_x / 90.0) * (7*y[i] + 32*y[i+1] + 12*y[i+2] +
													32*y[i+3] + 7*y[i+4]);
	}

	return integral;
}
