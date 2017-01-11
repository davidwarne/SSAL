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
 * @brief a probability mass function 
 * @detail uses the lookup method
 * @param n number of events
 * @param p scaled probabilities (e.g., P(X = k) = p[k]/p_sum)
 * @param p_sum sum of entries in p
 *
 * @return X sampled the given PMF
 */
unsigned int durngpmfs(int n,double *p,double p_sum)
{
    double u;
    double tmp,q;
    int k;
    u = DURAND;
    tmp = 0;
    q = p_sum*u;
    for (k=0;k<n;k++)
    {
        tmp += p[k];
        if (tmp > q)
        {
            break;
        }
    }
    return (k < n) ? k : n-1;
}
