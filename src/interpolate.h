/*
 * interpolate.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <math.h>
#include <stdlib.h>

float linearInterp(float* x, float* y, int length, float pt);
float splineInterp(float* x, float* y, int length, float pt);
int findUpperIndex(float* x, int upper_bound, int lower_bound, float pt);

#endif /* INTERPOLATE_H_ */
