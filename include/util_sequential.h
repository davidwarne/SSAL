/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2015  David J. Warne
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef  __UTIL_SEQUENTIAL_H_
#define  __UTIL_SEQUENTIAL_H_

#include <math.h>
#include <stdlib.h>

#define ONE_ON_RAND_MAX (1.0/((float)RAND_MAX))


/**
 * @brief structure to hold random number generator information
 * @details this structure is heavily affected by compile options
 * if the Intel Mathe kernel libraries (MKL) are use for example then
 * the Vector statistics library generator streams are used. 
 * The default however is to simply hold a seed and a function pointer
 * to a uniform RNG function (which is not considered thread safe).
 */
struct sRNG_struct {
#if defined(__MKL__)
#elif defined(__GSL__)
    /**seed value*/
    unsigned int seed;
    /**GSL random number generator*/
    gsl_rng * r;
    gsl_rng_type * T;
#else
    /**seed value*/
    unsigned int seed;
#endif
    /*seed function*/
    void (*s)(unsigned int);
    /** uniform(0,1) generator function*/
    int (*U)(void);
};
typedef struct sRNG_struct sRNG;

#define SURAND (((float)(*(__UTIL_sRNG.U))())/((float)RAND_MAX))
#define DURAND (((double)(*(__UTIL_sRNG.U))())/((double)RAND_MAX))

/** a global list of RNG streams */
extern sRNG __UTIL_sRNG;

void suarngs(unsigned int, void (*)(unsigned int), int (*)(void));

/* Continuous distributions samplers*/
float surngexps(float);
float surngus(float, float);

/*Discrete distribution samplers*/
unsigned int surngpmfs(int,float *,float);
unsigned int surngpois(float);

int suhzds(int ,int ,float * restrict, float *, float * restrict, float * restrict);

int sumlnls(int, int, float, float * restrict, float * restrict, float * restrict, 
    float * restrict, float, int, int, float, int, int * restrict, 
    int (*)(int float * float *), int * restrict);

#endif
