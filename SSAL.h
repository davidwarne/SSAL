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
 * @file SSAL.h
 * @brief A fast implementation of standard Stochastic Simulation algorithms.
 * @details This library provides optimised sequential and parallel exact and 
 * approximate stochastic simulation algorithms. The implementations provide an
 * easy to use application program interface (API) using standard C structures.
 *
 * @author David J. Warne (david.warne@qut.edu.au)
 * @author School of Mathematical Sciences
 * @author Science and Engineering Faculty
 * @author Queensland University of Technology
 *
 * @date 3 Sep 2015
 */

#ifndef __SSAL_H_
#define __SSAL_H_

#include <stdio.h>
#include <stdint.h>
#include "SSAL_internal.h"

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
    SSAL_ASSA_MULTI_LEVEL_SEQUENTIAL
};

/** type name for  SSAL_ChemicalReactionNetwork_struct*/
typedef struct SSAL_ChemicalReactionNetwork_struct SSAL_ChemicalReactionNetwork;

/** type name for  SSAL_RealisationSet_struct*/
typedef struct SSAL_RealisationSet_struct SSAL_RealisationSet;

/** type name for  SSAL_Algorithm_enum*/
typedef enum SSAL_Algorithm_enum SSAL_Algorithm_t;

/* function prototypes*/

/**
 * @brief Initilialises the SSAL library
 * @detail Initialised the computational backend based on compilation
 * options applied. Some user customisations are also available.
 *
 * @param argc the number of input args
 * @param argv an array of args
 *
 * @retVal SSAL_SUCCESS if initialisation was successful
 */
int SSAL_Inititalise(int,char **); 

/**
 * @brief Safely creates a chemical reaction network
 * @detail checks inputs are valid and creates the network
 *
 * @param N the number of chemical species
 * @param M the number of chemical reactions
 * @param nu_minus the reactant stochiometric coefficient matrix
 * @param nu_plus the product stochiometric coefficient matrix
 * @param c kinetic rate contstants
 *
 * @returns A Chemical Reaction Network Structure
 *
 * @note nu_minus and nu_minus must be N*M length vectors representing 
 * an M by N matrix using row-major storage. nu[j][i] is the quantity of chemical
 * species i required for reaction j.
 *
 * @note c must an M length vector, where c[j] is the kinetic rate of reaction j.
 *
 * @note Internally this function creates new copies of the stochiometric 
 * coefficient matrices. This ensures that the resulting allocated memory is valid, 
 * however it does not prevent invalid data in the matrices if the input pointers 
 * are not valid.
 */
SSAL_ChemicalReactionNetwork SSAL_CreateChemicalReactionNetwork(char ** names, int ,int, 
                                                        restrict float *, restrict  float *, float);
/**
 * @brief Performs N realisations of the a reaction network using the selected algorithm
 *
 * @param N Number of realisations to generate
 * @param model A chemical reaction network structure
 * @param X0 Initial chemical species copy numbers
 * @param T  A vector of time points to observe each realisation at
 * @param NT the number of time points 
 * @param alg the selected algorithm 
 * @param args algorithm specific parameters
 * 
 * @returns an array of realisations
 */
SSAL_RealisationSet SSAL_SimulateReactionNetwork(int, SSAL_ChemicalReactionNetwork, float *, float*, 
                                            int, SSAL_Algorithm, const char *);

/**
 * @briefs Estimates the expected values of copy numbers at give time points
 * 
 * @params model A chemical reaction network structure
 * @param X0 Initial chemical copy numbers
 * @param T  A vector of time points to observe each realisation at
 * @param NT the number of time points 
 * @param alg the selected algorithm 
 * @param args algorithm specific parameters
 */
SSAL_RealisationSet SSAL_ComputeExpectedValue(SSAL_ChemicalReactionNetwork, float *,float *, int,
                                        SSAL_Algorithm, const char *);

/**
 * @brief write simulation results to file
 * @param stream the output stream
 * @param data a Realisation set
 */
int SSAL_WriteRealisationSet(FILE *,SSAL_RealisationSet);


#endif
