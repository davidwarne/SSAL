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
/**
 * @file SSAL_internal.h
 * @brief Computational backend for SSAL
 * @detail Provides interfaces to al computational routines. All routines follow a 
 * standard interface convention which allows them all to be registered at runtime 
 * to a lookup table which is utilised by the API to access algorithms generically.
 * As a result switching algorithms is always just a matter of changing the flag.
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Science and Engineering Faculty
 * @author Queensland University of Technology
 *
 * @date 3 Sep 2015
 */

#ifndef __SSAL_INTERNAL_H_
#define __SSAL_INTERNAL_H_

#include "ESSA_sequential.h"
#include "ASSA_sequential.h"

#ifndef __PARALLEL__
#include "ESSA_parallel.h"
#include "ASSA_parallal.h"
#endif

/* Data structures and typedefs */

/**
 * @struct SSAL_ChemicalReactionNetwork_struct
 * @brief Defines a stochastic Chemical Reaction network model
 * @detail A user can manually interact with an instance of this structure
 * or use the SSAL_Create_ChemicalReactionNetwork function
 */
struct SSAL_ChemicalReactionNetwork_struct {
    /** Chemical Species Names*/
    char ** names;
    /** number of Chemical Species*/
    uint32_t N;
    /** number of reactions */
    uint32_t M;
    /** stochiometric coefficient matrix of the reactants */
    float *nu_minus;
    /** stochiometric coefficient matrix of the products */
    float *nu_plus;
    /** kinetic rate constants*/
    float *c;
};

/**
 * @struct SSAL_RealisationSet_struct
 * @brief defines a storage struct for simulation data
 * @detail stores multiple realisations of a model at specific 
 * points in time. This struction is also used to store expected
 * value data along with uncertainties. 
 */
struct SSAL_RealisationSet_struct{
    /** Chemical Species Names*/
    char ** names;
    /** number of Chemical Species*/
    uint32_t N;
    /** number of realisations*/
    uint32_t NR;
    /** number of time points*/
    uint32_t NT;
    /** time points*/
    float *T;
    /** The An NR*NT*N array containing copy numbers
     * @detail Represents a NR x NT x N matrix where
     * X_r[k][j][i] is the copy number of the i-th species 
     * at the j-th time step in the k-th realiseation. Matrix is row-major
     */
    float *X_r;
    /** The An NR*NT*N array containing confidence intervals
     * @detail Represents a NR x NT x N matrix where
     * X_r[k][j][i] is the confidence intervals of the i-th species 
     * at the j-th time step in the k-th realiseation. Matrix is row-major.
     * @note this is only utilised when each realisation is an expectation
     */
    float *X_conf;

}

/**
 * @enum SSAL_Algorithm_enum
 * @brief Enumeration of the available Exact and Approximate algorithms
 */
enum SSAL_Algorithm_enum {

    SSAL_ESSA_GILLESPIE_SEQUENTIAL,
    SSAL_ESSA_GIBSON_BRUCK_SEQUENTIAL,
    SSAL_ASSA_TAU_LEAP_SEQUENTIAL,
    SSAL_ASSA_TIME_DISCRETISE_SEQUENTIAL,
    SSAL_ASSA_MULTI_LEVEL_SEQUENTIAL,
#ifdef __PARALLEL__
#endif
    SSAL_NUM_ALG
};

/** type name for  SSAL_ChemicalReactionNetwork_struct*/
typedef struct SSAL_ChemicalReactionNetwork_struct SSAL_ChemicalReactionNetwork;

/** type name for  SSAL_RealisationSet_struct*/
typedef struct SSAL_RealisationSet_struct SSAL_RealisationSet;

/** type name for  SSAL_Algorithm_enum*/
typedef enum SSAL_Algorithm_enum SSAL_Algorithm_t;

/** put function registration stuff in here */

#endif
