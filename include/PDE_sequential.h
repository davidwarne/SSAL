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

#ifndef __PDE_SEQUENTIAL_H_
#define __PDE_SEQUENTIAL_H_
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* BTCS method for 1-d parabolic PDE with non-linear diffusion and sources*/
int
dbtcsfps(double, double, int, double* restrict, double * restrict, 
         double * restrict, void (*)(double *, int , double*, double *), 
         void (*)(double *, int, double *, double *), int, int, double, int, 
         double * restrict);
#endif /*__PDE_SEQUENTIAL_H_*/
