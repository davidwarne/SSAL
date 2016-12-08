/* SSAL: Stochastic Simulation Algorithm Library
 * Copyright (C) 2016  David J. Warne
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

#ifndef __ODE_SEQUENTIAL_H_
#define __ODE_SEQUENTIAL_H_


/* RK4 method*/
int srk4s(int, int, int, float * restrict, float * restrict, float * restrict, void (*f)(float *, unsigned int, float *, unsigned int, float, float *), int , int * restrict, float, float * restrict);
int drk4s(int, int, int, double * restrict, double * restrict, double * restrict, void (*f)(double *, unsigned int, double *, unsigned int, double, double *), int , int * restrict, double, double * restrict);
#endif
