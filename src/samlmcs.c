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
#include "ASSA_sequential.h"
#include "util_sequential.h"

/**
 * @brief Single precision Multi-level Monte Carlo method
 * @datails An efficient Monte Carlo estimator with low statistical bias against exact simulation.
 * This function computes for a tau leap approximation Z_L with tau = tau0*M^(-L),
 *      E[f(Z_L(T))] ~ Q_0 + sum_l=1^L Q_l + O(epsilon)
 * where
 *      Q_0 = 1/n0 sum_j=1^n0 f(Z_0^j(T)) with Z_0^j the j-th realisation of Z_0
 *      Q_l = 1/nl sum_j=1^nl f(Z_l^j(T)) - f(Z_{l-1}^j(T))
 *
 * The function is an unbiased esitmator for E[f(Z_L(T))] thus,
 *      E[f(X(T))] = E[f(Z_L(T))] + O(tauL)
 *
 * @param m number of reactions
 * @param n dimension of state vector
 * @param nt number of timesteps to sample at
 * @param T arr of length nt of sample times
 * @param X0 Initial condition
 * @param nu_minus reactant stoichiometric coefficients m x n (row-major format)
 * @param nu the stoichiometric matrix m x n (row-major format)
 * @param c the kinetic reaction rates
 * @param tau0 base level tau such that taul = tau0*M^(-l)
 * @param M scale factor between levels
 * @param L the number of levels
 * @param epsilon the desired statistical error (desired standard deviation)
 * @param ndims numebr dimension to measure
 * @param dims the dimensions indices 
 * @param E_X expected value of f(X(T))
 * @param V_X varience of f(X(T))
 * @param f functional of state vector f : \mathbb{Z}^n -> \mathbb{R}
 *
 * @param if f == NULL then the expected state vector is used
 *
 */
int samlmcbs(int m,int n, int nt, float * restrict T, float * restrict X0, float * restrict nu_minus,
    float * restrict nu, float * restrict c, float tau0, int M, int L, float epsilon, int ndims,int restrict dims, float * restrict E_X, float * restrict V_X, float (*f)(int,float *))
{

}
