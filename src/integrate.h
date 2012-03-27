/*
 * Integrator.h
 *
 *  Created on: Feb 9, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef INTEGRATE_H_
#define INTEGRATE_H_

#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>

/* Numerical integration methods */
typedef enum integrationSchemes {
	RIEMANN_LEFT,
	RIEMANN_RIGHT,
	RIEMANN_CENTER,
	TRAPEZOIDAL,
	SIMPSONS,
	SIMPSONS38,
	BOOLES,
	ENUM_END,
} integrationScheme;


float integrate(float* x, float* y, int length, integrationScheme scheme);
void cumulativeIntegral(float* x, float* y, float* cdf, int length,
												integrationScheme scheme);
float computeRiemannRight(float* x, float* y, int length);
float computeRiemannLeft(float* x, float* y, int length);
float computeRiemannCenter(float* x, float* y, int length);
float computeTrapezoidal(float* x, float* y, int length);
float computeSimpsons(float* x, float* y, int length);
float computeSimpsons38(float* x, float* y, int length);
float computeBooles(float* x, float* y, int length);

#endif /* INTEGRATE_H_ */
