
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

#ifndef __ASSA_SEQUENTIAL_H_
#define __ASSA_SEQUENTIAL_H_

#include "util_sequential.h"

#define MLMC_TIMING_TRIALS 100

/*Tau-leaping method*/
int 
datauls(int,int,int,double * restrict, double *restrict, double * restrict, 
        double * restrict, double *restrict, int, int * restrict, double,
        double * restrict);

/*Correlated tau-leaping method (for nested tau levels)*/
int 
dactauls(int, int, int, double * restrict, double * restrict, double * restrict, 
         double * restrict, double * restrict, int, int * restrict, double, int, 
         double * restrict, double * restrict);
/* Correlated tau-leaping method with Exact Next reaction method*/
int 
daectauls(int, int, int, double *restrict, double * restrict,double * restrict, 
          double *restrict, double * restrict, int, int *, double,
          double * restrict, double * restrict);


/* Euler-Maruyama Method*/
int 
daems(int , int, int , double * restrict , double * restrict, double * restrict, 
      void (*)(double *, unsigned int, double *, unsigned int, double,double *), 
      void (*)(double*, unsigned int, double *, unsigned int, double,double*),
      int , int *restrict , double , double *restrict );

/*Correlated Euler-Maruyama method (for nested h steps*/
int 
dacems(int,  int, int , double * restrict ,double * restrict, double * restrict, 
       void (*)(double*, unsigned int, double*, unsigned int, double, double*), 
       void (*)(double*, unsigned int, double *, unsigned int, double, double*), 
       int , int * restrict , double , int , double * , double * );

/* discrete-time lattice-based random walk with proliferation*/
int
dalrws(int, int *, void (*)(int, int *, int, int *,int *), int, int, 
       double * restrict, double * restrict, double * restrict, 
       double (*)(double, double *), double (*)(double, double*), 
       int, int * restrict, double, double *);
#endif /* __ASSA_SEQUENTIAL_H_ */
