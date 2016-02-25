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
 * @brief computes generic propensity function, or hazard function
 * for a discrete state continous time markov process
 * 
 *
 * @param m number of reactions
 * @param n dimension of state vector
 * @param R reactant stoichiometric coefficients m x n (stored in row-major format)
 * @param c rate constants
 * @param X state vector
 * @param A output propensities
 */
int duhzds(int m,int n,double * restrict R, double *restrict c, double * restrict X, double * restrict A)
{
    int i,j,k;
    /*compute 
     * a_j(X) = c_j \Pi_(i=1)^N nchoosek(X_i,R_j,i)*factorial(R_j,i) 
     * */
    for (j=0;j<m;j++)
    {
        A[j] = c[j];
        for(i=0;i<n;i++)
        {
            for (k=0;k<R[j*n + i];k++)
            {
                A[j] *= (X[i] - k);
            }
        }
    }
    return 0;
}

