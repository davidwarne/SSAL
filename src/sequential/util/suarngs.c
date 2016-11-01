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

sRNG __UTIL_sRNG;

#if defined(__MKL__)
int mkl_rand(void);
void mkl_srand(unsigned int);
#endif

/**
 * @brief adds a new RNG to the global 
 * @param seed a seed value
 * @param s function pointer to seed function
 * @param u function pointer to generator function
 *
 * @note if s and u are NULL then the default C functions srand and rand are used
 */
void suarngs(unsigned int seed, void (*s)(unsigned int), int (*u)(void))
{
#if defined(__MKL__)
    __UTIL_sRNG.seed = seed;
    __UTIL_sRNG.brng = VSL_BRNG_MT19937;
    __UTIL_sRNG.s = &mkl_srand;
    __UTIL_sRNG.U = &mkl_rand;
#elif defined(__GSL__)
   gsl_rng_env_setup();
   __UTIL_sRNG.T = gsl_rng_default;
   __UTIL_sRNG.r = gsl_rng_alloc(__UTIL_sRNG.T);
#else
    __UTIL_sRNG.seed = seed;
    if (s == NULL)
    {
        __UTIL_sRNG.s = &srand;
    }
    else
    {
        __UTIL_sRNG.s = s;
    }

    if (u == NULL)
    {
        __UTIL_sRNG.U = &rand;
    }
    else
    {
        __UTIL_sRNG.U = u;
    }
#endif
    (*(__UTIL_sRNG.s))(seed);
}

#if defined(__MKL__)
int mkl_rand(void)
{
    if (__UTIL_sRNG.ind < RNG_BLOCK_SIZE -1)
    {
        __UTIL_sRNG.ind++;
    }
    else
    {
        vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD,__UTIL_sRNG.stream,RNG_BLOCK_SIZE,__UTIL_sRNG.dbuf,0,1);
        __UTIL_sRNG.ind = 0;
    }
}
void mkl_srand(unsigned int s)
{
    vslNewStream(&(__UTIL_sRNG.stream),__UTIL_sRNG.brng,s);
    __UTIL_sRNG.ind = RNG_BLOCK_SIZE;
}
#endif
