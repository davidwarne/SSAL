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

#include "util_sequential.h"

/**
 * @define generates a exponential random variate
 * @define uses the inverse transform method
 * @param lambda the distribution parameter 
 *
 * @returns X ~ Exp(lambda)
 */
double durngexps(double lambda)
{
    int i;
#if defined(__MKL__)
#elif defined(__GSL__)
   return gsl_ran_exponential(__UTIL_sRNG.r,(double)lambda[i]);
#else
    double u;
    u = DURAND;
    return -log(u)/lambda;
#endif
}
