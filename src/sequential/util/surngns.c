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
 * @define generates a normal random variate
 * @define uses the Box-Muller method
 * @param mu the mean
 * @param sigma the standard deviation
 *
 * @returns X ~ N(mu,sigma)
 */
float surngns(float mu,float sigma)
{
    float u1,u2;
    u1 = SURAND;
    u2 = SURAND;
    return sigma*sqrtf(-2.0*logf(u1))*cosf(2.0*M_PI*u2);
}

