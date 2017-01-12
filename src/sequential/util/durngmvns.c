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

#include "util_sequential.h"

/**
 * @brief generate Mulit-variate normal random variates 
 * @details uses Box-Muller method for generating independent pairs
 * @param d dimension of distribution
 * @param mu mean vector
 * @param lambda lower triangular matrix (stored in row-major format) 
 *               such that lambda*lambda^T = Sigma where Sigma is the 
 *               covariance matrix.
 * @param Z N(0,I) random sample
 * @param X N(mu,Sigma) random sample
 */
void 
durngmvns(int d, double * restrict mu, double * restrict lambda, 
          double * restrict Z, double * restrict X)
{
    double u1,u2,r,t;
    int i,j;
    for (i=0;i<(d-1);i+=2)
    {
        /*generate indepentent standard normal pairs*/
        u1 = DURAND;
        u2 = DURAND;
        r = sqrt(-2.0*log(u1));
        t = 2.0*M_PI*u2;
        Z[i] = r*cos(t);
        Z[i+1] = r*sin(t);
    }
    /*check if d is odd, if so generate one more variate*/
    if (d%2 == 1)
    {
        u1 = DURAND;
        u2 = DURAND;
        r = sqrt(-2.0*log(u1));
        t = 2.0*M_PI*u2;
        Z[i] = r*cos(t);
    }
    
    /*compute X = lambda*Z + mu*/
    if (lambda != NULL)
    {
        for (i=0;i<d;i++)
        {
            X[i] = 0.0;
            for (j=0;j<d;j++)
            {
                X[i] += lambda[i*d +j]*Z[j];
            }
        }
    }
    else
    {
        for (i=0;i<d;i++)
        {
            X[i] = Z[i];
        }
    }

    if (mu != NULL)
    {
        for (i=0;i<d;i++)
        {
            X[i] += mu[i];
        }
    }
}
