
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

#ifndef __ASSA_SEQUENTIAL_H_
#define __ASSA_SEQUENTIAL_H_

#include "util_sequential.h"

#define MLMC_TIMING_TRIALS 100

/*Tau-leaping method*/
int satauls(int,int,int,float * restrict, float *restrict, float * restrict, 
    float * restrict, float *restrict, int, int * restrict, float, float * restrict);

int datauls(int,int,int,double * restrict, double *restrict, double * restrict, 
    double * restrict, double *restrict, int, int * restrict, double, double * restrict);

/*Correlated tau-leaping method (for nested tau levels)*/
int sactauls(int, int, int, float * restrict, float * restrict, float * restrict, 
    float * restrict, float * restrict, int, int * restrict, float, int, float * restrict, 
    float * restrict);

int dactauls(int, int, int, double * restrict, double * restrict, double * restrict, 
    double * restrict, double * restrict, int, int * restrict, double, int, double * restrict, 
    double * restrict);

/* Multilevel Monte Carlo estimator for DSCT Markov Processes (bias estimator)*/
int samlmcbs(int ,int , int , float * restrict , float * restrict , float * restrict ,
    float * restrict , float * restrict , float , int , int , float, 
    int,int * restrict , int (*)(int,float *, float *), float * restrict , float * restrict );

int damlmcbs(int ,int , int , double * restrict , double * restrict , double * restrict ,
    double * restrict , double * restrict , double , int , int , double, 
    int,int * restrict , int (*)(int,double *, double *), double * restrict , double * restrict );

/* Euler-Maruyama Method*/
int saems(int, int , int , float * restrict , float * restrict, float * restrict , void (*)(float *, unsigned int, float *, unsigned int , float,float *), void (*)(float*,unsigned int, float *, unsigned int, float,float*),int , int *restrict , float , float *restrict );

int daems(int , int, int , double * restrict , double * restrict, double * restrict , void (*)(double *, unsigned int, double *, unsigned int, double,double *), void (*)(double*, unsigned int, double *, unsigned int, double,double*),int , int *restrict , double , double *restrict );

/*Correlated Euler-Maruyama method (for nested h steps*/
int dacems(int,  int, int , double * restrict ,double * restrict, double * restrict , void (*)(double *, unsigned int, double *, unsigned int, double, double *), void (*)(double*, unsigned int, double *, unsigned int, double, double *), int , int * restrict , double , int , double * , double * );
#endif
