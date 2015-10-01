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

#include "util_sequential.h"

/**
 * @brief single precision Poisson random variates generator
 * @detail For small lambda Donald Knuth's method is used
 * otherwise the method proposed by Atkinson is used. 
 * If compiled with GSL or MKL then those libraries are utilised
 */
unsigned int surngpois(float lambda)
{
    int i;
#if defined(__MKL__)
#elif defined(__GSL__)
    return (float)gsl_ran_poisson(__UTIL_sRNG.r,(double)lambda[i]);
#else
   if (lambda < 30.0)
   {
       float L;
       int k;
       float p;
       L = expf(-lambda);
       k = 1;
       p = SURAND;
       while (p > L)
       {
           k++;
           p = p*SURAND;
       }
       return  k-1;
   }
   else
   {
       float c,k;
       float beta,alpha;
       c = 0.767 - 3.36/lambda;
       beta = M_PI/sqrtf(3.0*lambda);
       alpha = beta*lambda;
       k = logf(c) - lambda - logf(beta);

       while(1)
       {
           float u,v,x,y,lhs,rhs,t;
           unsigned int n;
           u = SURAND;
           x = (alpha - logf((1.0 - u)/u))/beta;
           n = floorf(x + 0.5);
           if (n < 0)
           {
               continue;
           }
           v = SURAND;
           y = alpha - beta*x;
           t = 1.0 + expf(y);
           lhs = y + logf(v/(t*t));
           rhs = k + n*logf(lambda) - lgammaf(n+1);
           if (lhs <= rhs)
           {
               return n;
           }
       }
   }
#endif
}
