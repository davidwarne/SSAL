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

#ifndef __SSAL_H_
#define __SSAL_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "SSAL_internal.h"
/*error codes*/
#define SSAL_SUCCESS            0
#define SSAL_INVALID_OPTION     -1
#define SSAL_UNKNOWN_OPTION     -2
#define SSAL_IO_ERROR           -3
#define SSAL_UNKNOWN_TYPE       -4

/*option codes*/
#define SSAL_ALG_SEQUENTIAL     0x01
#define SSAL_ALG_PARALLEL       0x02
#define SSAL_ALG_DISTRUBUTED    0x04
#define SSAL_ALG_ACCEL          0x08
#define SSAL_ALG_DEFAULT        0x0F

#define SSAL_MAX_NAME_SIZE      128

#define SSAL_REGISTER_ALG(alg,icf,ef,cpf,ovf) \
                        SSAL_registerInputCheckFunc((alg),(icf));\
                        SSAL_registerExecuteFunc((alg),(ef));\
                        SSAL_registerCheckpointFunc((alg),(cpf));\
                        SSAL_registerOutputValidationFunc((alg),(ovf));\
#define SSAL_DISABLE_ALG(alg) (SSAL_ALG_ENABLED[(alg)] = 0)  
#define SSAL_ENABLE_ALG(alg) (SSAL_ALG_ENABLED[(alg)] = 1)  

/**@struct SSAL_Options_struct
 * @brief Contains user runtime configurable options
 */
struct SSAL_Options_struct {
    /**bitmask of algorithm types that are available for selection*/
    uint8_t algOptions; 
};

typedef struct SSAL_Options_struct SSAL_Options;

/**
 * @enum SSAL_AlgorithmType_enum
 * @brief Enumeration of the available Exact and Approximate algorithms
 */
enum SSAL_AlgorithmType_enum {
    SSAL_ESSA_GILLESPIE_SEQUENTIAL,
    SSAL_ESSA_GIBSON_BRUCK_SEQUENTIAL,
    SSAL_ASSA_TAU_LEAP_SEQUENTIAL,
    SSAL_ASSA_TIME_DISCRETISE_SEQUENTIAL,
    SSAL_ASSA_MULTI_LEVEL_SEQUENTIAL,
    SSAL_ESSA_AUTO,
    SSAL_ASSA_AUTO,
    SSAL_RESTORE,
    SSAL_INFO,
    SSAL_USER_ALG,
    SSAL_NUM_ALGS
};
/** type name for  SSAL_Algorithm_enum*/
typedef enum SSAL_AlgorithmType_enum SSAL_AlgorithmType;

enum ModelType_enum {
    SSAL_CHEMICAL_REACTION_NETWORK
};
/** type name for  SSAL_ModelType_enum*/
typedef enum SSAL_ModelType_enum SSAL_ModelType;

enum SimulationType_enum {
    SSAL_REALISATIONS,
    SSAL_EXPECTEDVALUE
};
/** type name for  SSAL_SimulationType_enum*/
typedef enum SSAL_SimulationType_enum SSAL_SimulationType;

/**
 * generic stochastic model 
 */
struct SSAL_Model_struct {
    SSAL_ModelType type;
    void *model;
}

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
/** type name for  SSAL_ChemicalReactionNetwork_struct*/
typedef struct SSAL_ChemicalReactionNetwork_struct SSAL_ChemicalReactionNetwork;

/**
 * @struct SSAL_Simulation_struct
 * @brief defines generic simulation 
 */
struct SSAL_Simulation_struct{
    /**
     * The type of simulation structure
     */
    SSAL_SimulationType type;
    /** model to use in simulation */
    SSAL_model *model;
    void *sim;
};
/** type name for  SSAL_Simulation_struct*/
typedef struct SSAL_Simulation_struct SSAL_Simulation;

/**
 * realisation simulation
 */
struct SSAL_RealisationSimulation_struct {
    /**numbr of observation times */
    int NT;
    /** observation times */
    float *T;
    /** number of realisations*/
    int NR;
    /**number of initial conditions */
    int NIC;
    /** initial conditions */
    float *IC;
    /** number od vvariables to observe */
    int Nvar;
    /** model variables to observe by name */
    char **var;
    /** output data*/
    float *output;
};
/** type name for  SSAL_RealisationSimulation_struct*/
typedef struct SSAL_RealisationSimulation_struct SSAL_RealisationSimulation;


struct SSAL_Algorithm_struct {
    int (*inputCheck_func)(SSAL_Simulation *); 
    int (*execute_func)(SSAL_Simulation *,int, char **,int (*)(FILE *,void *),unsigned char,void *);
    int (*checkpoint_func)(FILE *,void *);
    int (*outputValidation_func)(Simulation *);    
};
typedef struct SSAL_Algorithm_struct SSAL_Algorithm;

#define SSAL_INPUTCHECK(alg,s) (*(SSAL_ALG_LUT[(alg)].inputCheck_func))((s))
#define SSAL_EXEC(alg,s,a,v) (*(SSAL_ALG_LUT[(alg)].execcute_func))((s),(a),(v),NULL,0,NULL)

/**
 * @brief Register a Algorithm Input Check Callback
 * @param alg the algorithm type
 * @param f input check function
 */
int SSAL_registerInputCheckFunc(SSAL_AlgorithmType,int (*)(SSAL_Simulation *));
/**
 * @brief Register a Algorithm Execute Callback
 * @param alg the algorithm type
 * @param f execution function
 */
int SSAL_registerExecuteFunc(SSAL_AlgorithmType,int, char **,
                        int (*)(SSAL_Simulation *,int (*)(FILE *,void *),unsigned char, void *));
/**
 * @brief Register a Algorithm Check Point Callback
 * @param alg the algorithm type
 * @param f check point function
 */
int SSAL_registerCheckPointFunc(SSAL_AlgorithmType,int (*)(FILE *,void *));
/**
 * @brief Register a Algorithm Output Validation Callback
 * @param alg the algorithm type
 * @param f output validation function
 */
int SSAL_registerOutputValidationFunc(SSAL_AlgorithmType,int (*)(FILE *,void *));


/* function prototypes*/

int SSAL_Inititalise(int,char **); 
SSAL_Model SSAL_CreateChemicalReactionNetwork(char **, int ,int, 
                            float * restrict nu_minus, float * restrict nu_plus, float * restrict c)
SSAL_Simulation SSAL_CreateRealisationsSim(SSAL_Model *,int, char **, int, int, 
                                float*, float * );
int SSAL_Simulate(SSAL_Simulation, SSAL_AlgorithmType, const char *);
int SSAL_WriteChemicalReactionNetwork(FILE *,SSAL_ChemicalReactionNetwork);
int SSAL_WriteSimulation(FILE *,SSAL_Simulation);
void SSAL_HandleError(int , char *,int, unsigned char, unsigned char, char *); 


#endif
