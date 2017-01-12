/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2017  David J. Warne
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

#if defined(__MKL__)
    #include <mkl.h>
    #include <mkl_vsl.h>
    #define RNG_BLOCK_SIZE 1000
#endif /* __MKL__ */

#define ONE_ON_RAND_MAX (1.0/((double)RAND_MAX))


/**
 * @brief structure to hold random number generator information
 * @details this structure is heavily affected by compile options
 * if the Intel Mathe kernel libraries (MKL) are use for example then
 * the Vector statistics library generator streams are used. 
 * The default however is to simply hold a seed and a function pointer
 * to a uniform RNG function (which is not considered thread safe).
 */
struct sRNG_struct 
{
#if defined(__MKL__)
    /** seed value*/
    unsigned int seed;
    /** MKL VSL random stream */
    VSLStreamStatePtr stream;
    /** Basic RNG type*/
    MKL_INT brng;
    /** Buffer of uniform Random variates
     * This is use as MKL functions ideally operate of arrays of length > 1000
     */
    double dbuf[RNG_BLOCK_SIZE];
    /** current position  in buffer */
    int ind; 
#elif defined(__GSL__)
    /**seed value*/
    unsigned int seed;
    /**GSL random number generator*/
    gsl_rng * r;
    gsl_rng_type * T;
#else /* not __MKL__ or __GSL__ */
    /**seed value*/
    unsigned int seed;
#endif /* not  __MKL__ or  __GSL__ */
    /**seed function*/
    void (*s)(unsigned int);
    /** uniform(0,1) generator function*/
    int (*U)(void);
};
typedef struct sRNG_struct sRNG;

#if defined(__MKL__)
    #define DURAND (__UTIL_sRNG.dbuf[(*(__UTIL_sRNG.U))()])
#else /* not __MKL__ */
    #define DURAND (((double)(*(__UTIL_sRNG.U))())/((double)RAND_MAX))
#endif /* not __MKL__ */

/** a global list of RNG streams */
extern sRNG __UTIL_sRNG;

void 
suarngs(unsigned int, void (*)(unsigned int), int (*)(void));

/* Continuous distributions samplers*/
double 
durngexps(double);

double 
durngus(double, double);

/* normal distribution*/
double 
durngns(double, double);

/* multi-variate normal distribution*/
void 
durngmvns(int , double * restrict, double * restrict , double * restrict , 
          double * restrict );

/*Discrete distribution samplers*/
unsigned int 
durngpmfs(int,double *,double);

unsigned int 
durngpois(double);

/*hazard rate/ propensity function computation*/
int 
duhzds(int ,int ,double * restrict, double *, double * restrict, 
       double * restrict);
#endif
