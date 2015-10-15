
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

int satauls(int,int,int,float * restrict, float *restrict, float * restrict, 
    float * restrict, float *restrict, int, int * restrict, float, float * restrict);

int sactauls(int, int, int, float * restrict, float * restrict, float * restrict, 
    float * restrict, float * restrict, int, int * restrict, float, int, float * restrict, 
    float * restrict);


int samlmcbs(int ,int , int , float * restrict , float * restrict , float * restrict ,
    float * restrict , float * restrict , float , int , int , float, 
    int,int * restrict , int (*)(int,float *, float *), float * restrict , float * restrict );

#endif
