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
#include <stdio.h>
/**
 * @brief single precision Poisson random variates generator
 * @detail For small lambda Donald Knuth's method is used
 * otherwise the method proposed by Atkinson is used. 
 */
unsigned int durngpois(double lambda)
{
    int i;
    if (lambda < 30.0)
    {
       double L;
       int k;
       double p;
       L = exp(-lambda);
       k = 1;
       p = DURAND;
       while (p > L)
       {
           k++;
           p = p*DURAND;
       }
       return  k-1;
   }
   else
   {
       double c,k;
       double beta,alpha;
       c = 0.767 - 3.36/lambda;
       beta = M_PI/sqrt(3.0*lambda);
       alpha = beta*lambda;
       k = log(c) - lambda - log(beta);

       while(1)
       {
           double u,v,x,y,lhs,rhs,t;
           unsigned int n;
           u = DURAND;
           x = (alpha - log((1.0 - u)/u))/beta;
           n = floor(x + 0.5);
           if (n < 0)
           {
               continue;
           }
           v = DURAND;
           y = alpha - beta*x;
           t = 1.0 + exp(y);
           lhs = y + log(v/(t*t));
           rhs = k + n*log(lambda) - lgammaf(n+1);
           if (lhs <= rhs)
           {
               return n;
           }
       }
   }
}
