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
 * @file SSAL.c
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
 * @date 20 Sep 2015
 */

#include "SSAL.h"

/**
 * @brief Handle error reporting
 * @param errCode The Error code 
 * @param funcName The Function namee in which the error occurred
 * @param lineNum The line number at which the error occurred
 * @param fatal Flag indicating the error was fatal and the progam should be aborted
 * @param warning Indicates the error is to be interpreted as a warning
 * @param msg A custom message which can be appended to the default
 */
void SSAL_HandleError(int errCode, char * funcName, int lineNum, unsigned char fatal, 
                    unsigned char warning, char *msg)
{
    char * errName;
    char * defaultMsg;
    char * errType;
    switch (errCode)
    {
        case SSAL_UNKNOWN_OPTION_ERROR:
            errName = "SSAL_UNKNOWN_OPTION_ERROR";
            defaultMsg = "A configuration option was not recognised.";
            break;
        default:
            errName = "SSAL_UNKNOWN_ERROR";
            defaultMsg = "No idea what happened then...";
            break;
    }

    if (warning)
    {
        errType = "ERROR";   
    }
    else
    {
        errType = "WARNING";
    }
        
    if (msg != NULL)
    {
        fprintf(stderr,"[%s %s]: %s\n %s\n Error Code [%d] Line [%d]\n",
            errType,errName,defaultMsg,msg,errCode,lineNum);    
    }
    else
    {
        fprintf(stderr,"[%s %s]: %s\n Error Code [%d] Line [%d]\n",
            errType,errName,defaultMsg,errCode,lineNum);    
    }
    

    if (fatal)
    {
        fprintf(stderr,"Fatal Error Aborting!\n");
        exit(1);
    }
    return;
}



/**
 * @brief Initilialises the SSAL library
 * @detail Initialised the computational backend based on compilation
 * options applied. Some user customisations are also available.
 *
 * @param argc the number of input args
 * @param argv an array of args
 * @todo Will add more options here as needed.
 * @retVal SSAL_SUCCESS if initialisation was successful
 */
int SSAL_Initialise(int argc,char **argv)
{
    int i;
    unsigned int seed;
    seed = 1337;
    /*search for config options*/
    for (i=1;i<argc;i++)
    {
        if (!strcmp("--SSALOPT",argv[i]))
        {
            i++;
            break;
        }
    }
    

    /*set defaults*/
    /*parse args*/
    for (;i<argc;i++)
    {
        if (!strcmp("--seed",argv[i]))
        {
            seed = (unsigned int)atoi(argv[++i]);          
        }
        else
        {
            SSAL_HandleError(SSAL_UNKNOWN_OPTION_ERROR,"SSAL_Initialise",__LINE__,0,1,argv[i]);
        }
    }

    /*initialise RNGS */
    #ifdef __SERIAL__
       suarngs(seed,NULL,NULL);
    #endif
}



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
SSAL_Model SSAL_CreateChemicalReactionNetwork(char ** names, int M,int N, 
                            float * restrict nu_minus, float * restrict nu_plus, float * restrict c)
{
    int i;
    char * funcname;
    char * nameArray;
    SSAL_Model newModel;
    SSAL_ChemicalReactionNetwork *newCRN;
   
   funcname = "SSAL_CreateChemicalReactionNetwork";
   
    /*build wrapper struct*/
    newModel.type = SSAL_CHEMICAL_REACTION_NETWORK;
    if ((newCRN = (SSAL_ChemicalReactionNetwork*)malloc(sizeof(SSAL_ChemicalReactionNetwork)))== NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
   
    newModel.model = (void *)newCRN;
    
    /*allocate memory*/
    newCRN->N = (uint32_t)N;
    newCRN->M = (uint32_t)M;
    
    if ((nameArray = (char *)malloc(newCRN->N*SSAL_MAX_NAME_SIZE*sizeof(char))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((newCRN->names = (char **)malloc(newCRN->N*sizeof(char*)))== NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    
    for (i=0;i<newCRN->N;i++)
    {
        newCRN->names[i] = nameArray + i*SSAL_MAX_NAME_SIZE;
    }
    
    if((newCRN->nu_minus = (float *)malloc(newCRN->N*newCRN->M*sizeof(float))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((newCRN->nu_plus = (float *)malloc(newCRN->N*newCRN->M*sizeof(float))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((newCRN->c = (float *)malloc(newCRN->M*sizeof(float))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    
    /*copy data to valid arrays*/
    for (i=0;i<newCRN->N*newCRN->M;i++)
    {
        newCRN->nu_minus[i] = nu_minus[i];
    }
    for (i=0;i<newCRN->N*newCRN->M;i++)
    {
        newCRN->nu_plus[i] = nu_plus[i];
    }
    for (i=0;i<newCRN->M;i++)
    {
        newCRN->c[i] = c[i];
    }
    for (i=0;i<newCRN->N;i++)
    {
        int j;
        for (j=0;j<SSAL_MAX_NAME_SIZE-1;j++)
        {
            newCRN->names[i][j] = names[i][j];
            if (names[i][j] == '\0')
            {
                continue;
            }
        }
        newCRN->names[i][SSAL_MAX_NAME_SIZE-1] = '\0';
    }

    /*all done, return the new model*/
    return newModel;
}

/**
 * @brief write chemical reation network too file
 * @param stream the output stream
 * @param data a chemical reaction network
 */
int SSAL_WriteChemicalReactionNetwork(FILE * stream ,SSAL_ChemicalReactionNetwork data)
{
    int i,j;

    if (stream == NULL)
    {
        return SSAL_IO_ERROR;
    }

    /*print species data*/
    fprintf(stream,"Chemical Species:\n");
    for (i=0;i<data.N;i++)
    {
        fprintf(stream,"\tX[%d]: %s\n",i,data.names[i]);
    }

    /*print equations*/
    fprintf(stream,"Chemical Reactions:\n");
    fprintf(stream,"\t         c \n");
    fprintf(stream,"\t nu^-*X --> nu^+*X\n");
    fprintf(stream,"\tnu^- =\n");
    for (j=0;j<data.M;j++)
    {
        fprintf(stream,"\t\t |");
        for (i=0;i<data.N;i++)
        {
            fprintf(stream," %f",data.nu_minus[j*data.N + i]);
        }
        fprintf(stream," |\n");
    }
    fprintf(stream,"\tnu^+ =\n");
    for (j=0;j<data.M;j++)
    {
        fprintf(stream,"\t\t |");
        for (i=0;i<data.N;i++)
        {
            fprintf(stream," %f",data.nu_plus[j*data.N + i]);
        }
        fprintf(stream," |\n");
    }
    fprintf(stream,"Kinetic Reaction Rates:\n");
    for (i=0;i<data.M;i++)
    {
        fprintf(stream,"\tc[%d] = %f\n",i,data.c[i]);
    }
    return SSAL_SUCCESS; 
}

/**
 * @brief write the simulation results to file
 * @param stream the output stream
 * @param sim a Simulation struct
 */
int SSAL_WriteSimulation(FILE *stream, SSAL_Simulation sim)
{
    int rc;
    char * funcname;
    funcname = "SSAL_WriteSimulation";
    switch (sim.type)
    {
        case SSAL_REALISATIONS:
            rc = SSAL_WriteRealisationsSim(stream,sim.sim);
            break;
        case SSAL_EXPECTEDVALUE:
            rc = SSAL_WriteExpectedValueSim(stream,sim.sim);
            break;
        default:
            rc = SSAL_UNKNOWN_TYPE_ERROR;
            break;
    }

    if (rc < 0)
    {
        SSAL_HandleError(rc,funcname,__LINE__,0,1,NULL);
    }
}

/**
 * @brief Write expected value data to a file
 * @param stream the output stream
 * @param sim the expected value simulation to export
 */
int SSAL_WriteExpectedValueSim(FILE * stream, SSAL_ExpectedValueSimulation * sim)
{
    int i,j,r;
    if (stream == NULL)
    {
        return SSAL_IO_ERROR;
    }

    fprintf(stream,"\"Time\"");
    for (i=0;i<sim->NT;i++)
    {
        fprintf(stream,",%f",sim->T[i]);
    }
    fprintf(stream,"\n");

    for (j=0;j<sim->Nvar;j++)
    {
        fprintf(stream,"\"E[%s]\"",sim->var[j]);
        for (i=0;i<sim->NT;i++)
        {
            fprintf(stream,",%f",sim->E[j*sim->NT + i]);
        }
        fprintf(stream,"\n");
        
        fprintf(stream,"\"V[%s]\"",sim->var[j]);
        for (i=0;i<sim->NT;i++)
        {
            fprintf(stream,",%f",sim->V[j*sim->NT + i]);
        }
        fprintf(stream,"\n");
    }
}

/**
 * @brief Write realisation data to a file
 * @param stream the output stream
 * @param sim the realisation sim to export
 */
int SSAL_WriteRealisationsSim(FILE * stream, SSAL_RealisationSimulation * sim)
{
    int i,j,r;
    if (stream == NULL)
    {
        return SSAL_IO_ERROR;
    }

    fprintf(stream,"\"Time\"");
    for (i = 0;i<sim->NT;i++)
    {
       fprintf(stream,",%f",sim->T[i]);
    }
    fprintf(stream,"\n");

    for (r=0;r<sim->NR;r++)
    {
       for (j=0;j<sim->Nvar;j++)
        {
            fprintf(stream,"\"%s\"",sim->var[j]);
            for (i=0;i<sim->NT;i++)
            {
                fprintf(stream,",%f",sim->output[r*(sim->Nvar*sim->NT) + j*(sim->NT) + i]);
            }
            fprintf(stream,"\n");
        }
    }
}

/**
 * @brief Creates a realisation simulation structure
 * @detail pre-allocates memory in an efficient structure for the model and algorithm
 * @param model a pointer referencing the model realisations will be made of
 * @param N number of variables to observe
 * @param obs a pointer of model entities to observe
 * @param NR number of realisations
 * @param NT number of time points to record
 * @param T time points to record
 * @param initCond initial conditions 
 */
SSAL_Simulation SSAL_CreateRealisationsSim(SSAL_Model *model, int N,char **obs, int NR, int NT, 
                                float* T, float * initCond)
{
    int i;
    char * obsArray;
    SSAL_Simulation newSim;
    SSAL_RealisationSimulation *newRS;
    char * funcname;
    funcname = "SSAL_CreateRealisationsSim";

    newSim.type = SSAL_REALISATIONS;
    newSim.model = model;
    
    if((newRS = (SSAL_RealisationSimulation*)malloc(sizeof(SSAL_RealisationSimulation)))==NULL)
    {
        SSAL_HandleError(SSAL_IO_ERROR,funcname,__LINE__,1,0,NULL);
    }

    newSim.sim = (void*)newRS;

    /*allocate memory*/
    newRS->NR = NR;
    newRS->NT = NT;
    
    if((newRS->T = (float *)malloc(NT*sizeof(float)))==NULL)
    {
        SSAL_HandleError(SSAL_IO_ERROR,funcname,__LINE__,1,0,NULL);
    }
    /**@todo extend to handle multiple scenarios */
    newRS->NIC = 1; 
    newRS->Nvar = N;
    if((obsArray = (char *)malloc(N*SSAL_MAX_NAME_SIZE*sizeof(char)))==NULL)
    {
        SSAL_HandleError(SSAL_IO_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((newRS->var = (char **)malloc(N*sizeof(char*)))==NULL)
    {
        SSAL_HandleError(SSAL_IO_ERROR,funcname,__LINE__,1,0,NULL);
    }
    for (i=0;i<newRS->Nvar;i++)
    {
        newRS->var[i] = obsArray + i*SSAL_MAX_NAME_SIZE;
    }

    /*data is always un-allocated first as it is model dependent*/
    newRS->output = NULL;

    /*copy data */
    for (i=0;i<newRS->NT;i++)
    {
        newRS->T[i] = T[i];
    }
    for (i=0;i<newRS->Nvar;i++)
    {
        int j;
        for (j=0;j<SSAL_MAX_NAME_SIZE-1;j++)
        {
            newRS->var[i][j] = obs[i][j];
            if (obs[i][j] == '\0');
            {
                continue;
            }
        }
        newRS->var[i][SSAL_MAX_NAME_SIZE-1] = '\0';
    }


    /*setting initial condition sizes are dependent on the model*/
    switch (model->type)
    {
        case SSAL_CHEMICAL_REACTION_NETWORK:
        {
            SSAL_ChemicalReactionNetwork *CRN;
            CRN = (SSAL_ChemicalReactionNetwork *)model->model;
            
            /*initial conditions will be the initial chemical species copy numbers*/
            newRS->IC = (float *)malloc((newRS->Nvar)*(CRN->N)*sizeof(float));
            for (i=0; i<(newRS->Nvar)*(CRN->N); i++)
            {
                newRS->IC[i] = initCond[i];
            }

            newRS->output = (float *)malloc((newRS->Nvar)*(newRS->NT)*(newRS->NR)*sizeof(float));
        }
            break;
        default:
            SSAL_HandleError(SSAL_UNKNOWN_TYPE_ERROR,funcname, __LINE__,1,0,NULL);        
            break;
    }
    
    return newSim;
}

/**
 * @brief Creates an Expected Value simulation structure
 * @detail pre-allocates memory in an efficient structure for the model and algorithm
 * @param model a pointer referencing the model realisations will be made of
 * @param N number of variables to observe
 * @param obs a pointer of model entities to observe
 * @param NR number of realisations
 * @param NT number of time points to record
 * @param T time points to record
 * @param initCond initial conditions 
 */
SSAL_Simulation SSAL_CreateExpectedValueSim(SSAL_Model *model, int N,char **obs, int NR, int NT, 
                                float* T, float * initCond)
{
    int i;
    char * obsArray;
    SSAL_Simulation newSim;
    SSAL_ExpectedValueSimulation *newEVS;
    newSim.type = SSAL_EXPECTEDVALUE;
    newSim.model = model;
    newEVS = (SSAL_ExpectedValueSimulation*)malloc(sizeof(SSAL_ExpectedValueSimulation));
    newSim.sim = (void*)newEVS;
    /*allocate memory*/
    newEVS->NR = NR;
    newEVS->NT = NT;
    newEVS->T = (float *)malloc(NT*sizeof(float));

    /**@todo extend to handle multiple scenarios */
    newEVS->NIC = 1; 
    newEVS->Nvar = N;
    obsArray = (char *)malloc(N*SSAL_MAX_NAME_SIZE*sizeof(char));
    newEVS->var = (char **)malloc(N*sizeof(char*));
    for (i=0;i<newEVS->Nvar;i++)
    {
        newEVS->var[i] = obsArray + i*SSAL_MAX_NAME_SIZE;
    }

    /*data is always un-allocated first as it is model dependent*/
    newEVS->E = NULL;
    newEVS->V = NULL;
    /*copy data */
    for (i=0;i<newEVS->NT;i++)
    {
        newEVS->T[i] = T[i];
    }
    for (i=0;i<newEVS->Nvar;i++)
    {
        int j;
        for (j=0;j<SSAL_MAX_NAME_SIZE-1;j++)
        {
            newEVS->var[i][j] = obs[i][j];
            if (obs[i][j] == '\0');
            {
                continue;
            }
        }
        newEVS->var[i][SSAL_MAX_NAME_SIZE-1] = '\0';
    }


    /*setting initial condition sizes are dependent on the model*/
    switch (model->type)
    {
        case SSAL_CHEMICAL_REACTION_NETWORK:
        {
            SSAL_ChemicalReactionNetwork *CRN;
            CRN = (SSAL_ChemicalReactionNetwork *)model->model;
            
            /*initial conditions will be the initial chemical species copy numbers*/
            newEVS->IC = (float *)malloc((newEVS->Nvar)*(CRN->N)*sizeof(float));
            for (i=0; i<(newEVS->Nvar)*(CRN->N); i++)
            {
                newEVS->IC[i] = initCond[i];
            }

            newEVS->E = (float *)malloc((newEVS->Nvar)*(newEVS->NT)*sizeof(float));
            newEVS->V = (float *)malloc((newEVS->Nvar)*(newEVS->NT)*sizeof(float));
        }
            break;
        default:
            SSAL_HandleError(SSAL_UNKNOWN_TYPE_ERROR,"SSAL_CreateRealisationsSimulation",
                            __LINE__,1,0,NULL);        
            break;
    }
    
    return newSim;
}


/**
 * @brief Runs a stochastic simulation
 *
 * @param sim pointer to simulation structure
 * @param alg the selected algorithm 
 * @param args algorithm specific parameters
 * 
 * @return return code
 */
int SSAL_Simulate(SSAL_Simulation *sim, SSAL_AlgorithmType alg, const char * args)
{
    int rc;
    switch(sim->model->type)
    {   
        case SSAL_CHEMICAL_REACTION_NETWORK:
            rc = SSAL_SimulateCRN(sim,alg,args);
            break;
        default:
            break;
    }
}

/**
 * @brief Run a realisation simulation
 * @param sim a Realisation simulation structure
 * @param alg the selected algorithm type
 * @param args algorithm specific args
 */
int SSAL_SimulateCRN(SSAL_Simulation *sim, SSAL_AlgorithmType alg, 
                        const char *args)
{
    int argc;
    char **argv;
    int rc;
    if (args == NULL)
    {
        argc = 0;
        argv = NULL;
    }
    else
    {
        argv = SSAL_UtilTokeniseArgs(&argc,args);
    }
    switch (sim->type)
    {
        case SSAL_REALISATIONS:
            rc = SSAL_SimulateCRNRealisations(sim->sim,sim->model->model,alg,argc,argv);
            break;
        case SSAL_EXPECTEDVALUE:
            rc = SSAL_SimulateCRNExpectedValue(sim->sim,sim->model->model,alg,argc,argv);
            break;
        default:
            break;
    }
    /*the first pointer is a pointer to the whole array*/
    if (args != NULL)
    {
        free(argv[0]);
    }
    return rc;
}

char * SSAL_GetArg(char *key,int argc,char **argv)
{
    int i;
    for (i=0;i<argc;i++)
    {
        if (!strcmp(key,argv[i]))
        {
            return argv[++i];
        }
    }
    return "";
}

/**
 * @brief run a realisation simulation on a Chemical Reaction Network model
 * @param sim a Realisation Simulation structure
 * @param model a Chemical Reaction network model
 * @param alg the selected algorithm type
 * @param argc the number of args
 * @param argv input ags
 */
int SSAL_SimulateCRNRealisations(SSAL_RealisationSimulation *sim, 
            SSAL_ChemicalReactionNetwork *model, SSAL_AlgorithmType alg, int argc, char ** argv)
{
    int j,i;
    int nvar;
    int var[sim->Nvar];
    float nu[model->N*model->M];
    float * X_rj;
    int varInd[sim->Nvar];

    for (i=0;i<sim->Nvar;i++)
    {
        for (j=0;j<model->N;j++)
        {
            if (!strcmp(sim->var[i],model->names[j]))
            {
                varInd[i] = j;
            }
        }
    }

    for (i=0;i<model->N*model->M;i++)
    {
        nu[i] = model->nu_plus[i] - model->nu_minus[i] ;
    }

    X_rj = sim->output;

    /*algorithm selector*/ 
    switch(alg)
    {
        case SSAL_ESSA_GILLESPIE_SEQUENTIAL:
        {
            for (j=0;j<sim->NR;j++)
            {
                segils(model->M,model->N,sim->NT,sim->T,sim->IC,model->nu_minus,
                    nu,model->c,sim->Nvar,varInd,X_rj+j*(sim->NT*sim->Nvar));
            }
        }
            break;
        case SSAL_ASSA_TAU_LEAP_SEQUENTIAL:
        {
            float tau;
            tau = (float)atof(SSAL_GetArg("--tau",argc,argv));

            for (j=0;j<sim->NR;j++)
            {
                satauls(model->M,model->N,sim->NT,sim->T,sim->IC,model->nu_minus,
                    nu,model->c,sim->Nvar,varInd,tau,X_rj+j*(sim->NT*sim->Nvar));
            }
        }
            break;
        default:
            break;
    }

}
 
 /**
 * @brief Compute Expected value on a Chemical Reaction Network model
 * @param sim an Expected Value Simulation structure
 * @param model a Chemical Reaction network model
 * @param alg the selected algorithm type
 * @param argc the number of args
 * @param argv input ags
 */
int SSAL_SimulateCRNExpectedValue(SSAL_ExpectedValueSimulation *sim, 
            SSAL_ChemicalReactionNetwork *model, SSAL_AlgorithmType alg, int argc, char ** argv)
{
    int j,i;
    int nvar;
    int var[sim->Nvar];
    float nu[model->N*model->M];
    float * X_r;
    double * E_X;
    double * V_X;
    int varInd[sim->Nvar];

    for (i=0;i<sim->Nvar;i++)
    {
        for (j=0;j<model->N;j++)
        {
            if (!strcmp(sim->var[i],model->names[j]))
            {
                varInd[i] = j;
            }
        }
    }

    for (i=0;i<model->N*model->M;i++)
    {
        nu[i] = model->nu_plus[i] - model->nu_minus[i] ;
    }

    X_r = (float *)malloc(sim->NT*sim->Nvar*sizeof(float));
    E_X = (double *)malloc(sim->NT*sim->Nvar*sizeof(double));
    V_X = (double *)malloc(sim->NT*sim->Nvar*sizeof(double));
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        E_X[i] = 0;
    }
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        E_X[i] = 0;
    }
    /*algorithm selector*/ 
    switch(alg)
    {
        case SSAL_ESSA_GILLESPIE_SEQUENTIAL:
        {
            /**@note Currently the Expected value computation is not ideal numerically
             * when a large number of realisations are used*/                
            for (j=0;j<sim->NR;j++)
            {
                segils(model->M,model->N,sim->NT,sim->T,sim->IC,model->nu_minus,
                    nu,model->c,sim->Nvar,varInd,X_r);
                
                for (i=0;i<sim->Nvar*sim->NT;i++)
                {
                    E_X[i] += X_r[i];
                }
                for (i=0;i<sim->Nvar*sim->NT;i++)
                {
                    V_X[i] += X_r[i]*X_r[i];
                }
            }
        }
            break;
        case SSAL_ASSA_TAU_LEAP_SEQUENTIAL:
        {
            float tau;
            tau = (float)atof(SSAL_GetArg("--tau",argc,argv));
            
            for (j=0;j<sim->NR;j++)
            {
                satauls(model->M,model->N,sim->NT,sim->T,sim->IC,model->nu_minus,
                    nu,model->c,sim->Nvar,varInd,tau,X_r);
                for (i=0;i<sim->Nvar*sim->NT;i++)
                {
                    E_X[i] += (double)X_r[i];
                }
                for (i=0;i<sim->Nvar*sim->NT;i++)
                {
                    V_X[i] +=(double)( X_r[i]*X_r[i]);
                }
            }
        }
            break;
        case SSAL_ASSA_MULTI_LEVEL_SEQUENTIAL:
        {
        }
            break;
        default:
            break;
    }

    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        E_X[i] /= (double)(sim->NR);
    }
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        V_X[i] /= (double)(sim->NR);
    }
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        V_X[i] -= E_X[i]*E_X[i];
    }

    /*convert variance in X to variance in the mean esitmator
     * i.e., Var[E[X(T)]] = Var[X]/n
     */
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        V_X[i] /= (double)(sim->NR);
    }

    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        sim->E[i] =  (float)E_X[i];
        sim->V[i] =  (float)V_X[i];
    }
}
 

/**
 * @brief utilitiy function which breaks a char array into a array of char arrays
 * @param argc arg count
 * @param args input arg string 
 * @return argv a array of pointers to a char arrays
 */
char** SSAL_UtilTokeniseArgs(int *argc,const char * args)
{
    char *buf;
    char **argv;

    int numChars;
    int numArgs;
    int argsLength;
    int i,j,k;
    char cur_char;
    char prev_char;

    cur_char = args[0];
    prev_char = ' ';
    numChars = 0;
    numArgs = 0;
    i = 0;

    /*first pass to get sizes... a bit insecure as buffer overruns can occur*/
    while (cur_char != '\0' && i != SSAL_MAX_BUFFER_SIZE)
    {
        numChars += (cur_char != ' ');
        numArgs += (cur_char != ' ' && prev_char == ' ');
        i++;
        prev_char = cur_char;
        cur_char = args[i];
    }


    argsLength = i; 
    buf = (char *)malloc((numChars+numArgs)*sizeof(char));
    argv = (char **)malloc(numArgs*sizeof(char*));
    j = 0;
    k = 0;
    /*second pass to collect tokenised data*/
    for (i=0;i<argsLength;i++)
    {
        /*starting a token*/
        if (args[i] == ' ' && args[i+1] != ' ')
        {
            j++;
            argv[k] = buf + j;
            k++;
        } /*ending a token*/
        else if (i==0 && args[i] != ' ')
        {
            buf[j] = args[i];
            argv[k] = buf + j;
            j++;
            k++;
        }
        else if ( args[i] != ' ')
        {
            buf[j] = args[i];
            j++;
        }
    }

    /*append null characters*/
    for (k=1;k<numArgs;k++)
    {
        argv[k][-1] = '\0';
    }
    buf[numChars+numArgs - 1] = '\0';
    *argc = numArgs;
    return argv;
}

