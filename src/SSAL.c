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
 * @date 12 Feb 2016
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
        case SSAL_UNSUPPORTED_ALGORITHM_ERROR:
            errName = "SSAL_UNSUPPORTED_ALGORITHM_ERROR";
            defaultMsg = "The selected simulation algorithm is not supported for this model.";
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
    unsigned char infoflag;
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
        else if (!strcmp("--info",argv[i]))
        {
            infoflag = 1;
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
    if (infoflag)
    {
        fprintf(stderr,"RNG Seed: %d\n",seed);
        fprintf(stderr,"RNG RAND_MAX: %d\n",RAND_MAX);
    }
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
                            SSAL_real_t * restrict nu_minus, SSAL_real_t * restrict nu_plus, SSAL_real_t * restrict c)
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
    
    if((newCRN->nu_minus = (SSAL_real_t *)malloc(newCRN->N*newCRN->M*sizeof(SSAL_real_t))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((newCRN->nu_plus = (SSAL_real_t *)malloc(newCRN->N*newCRN->M*sizeof(SSAL_real_t))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((newCRN->c = (SSAL_real_t *)malloc(newCRN->M*sizeof(SSAL_real_t))) == NULL)
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
        }
        newCRN->names[i][SSAL_MAX_NAME_SIZE-1] = '\0';
    }

    /*all done, return the new model*/
    return newModel;
}

/**
 * @brief Safely creates a multivariate stochastic differential equation
 * @detail checks inputs are valid and creates the SDE
 *
 * @param M the number of parameters  
 * @param N  the number of equations  
 * @param mu function pointer to drift function
 * @param sigma function pointer diffusion function
 * @param p model parameter vector
 *
 * @returns A Stochasti Differential Equation Structure
 */
SSAL_Model SSAL_CreateStochasticDifferentialEquation(char ** names, int M, int N, 
                            void (*mu)(SSAL_real_t *, uint32_t, SSAL_real_t*, uint32_t, SSAL_real_t,SSAL_real_t*), void (*sigma)(SSAL_real_t *, uint32_t, SSAL_real_t*, uint32_t, SSAL_real_t,SSAL_real_t*), SSAL_real_t *p)
{
    int i;
    char * funcname;
    char * nameArray;
    SSAL_Model newModel;
    SSAL_StochasticDifferentialEquation *newSDE;
   
   funcname = "SSAL_CreateStochasticDifferentialEquation";
   
    /*build wrapper struct*/
    newModel.type = SSAL_STOCHASTIC_DIFFERENTIAL_EQUATION;
    if ((newSDE = (SSAL_StochasticDifferentialEquation*)malloc(sizeof(SSAL_StochasticDifferentialEquation))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
   
    newModel.model = (void *)newSDE;
    
    /*allocate memory*/
    newSDE->N = (uint32_t)N;
    newSDE->M = (uint32_t)M;
    
    if ((nameArray = (char *)malloc(newSDE->N*SSAL_MAX_NAME_SIZE*sizeof(char))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    if((newSDE->names = (char **)malloc(newSDE->N*sizeof(char*)))== NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }
    
    for (i=0;i<newSDE->N;i++)
    {
        newSDE->names[i] = nameArray + i*SSAL_MAX_NAME_SIZE;
    }

    if ((newSDE->p = (SSAL_real_t *)malloc(newSDE->M*sizeof(SSAL_real_t))) == NULL)
    {
        SSAL_HandleError(SSAL_MEMORY_ERROR,funcname,__LINE__,1,0,NULL);
    }

    for (i=0;i<newSDE->M;i++)
    {
       newSDE->p[i] = p[i]; 
    }
   
    newSDE->mu = mu;
    newSDE->sigma = sigma;
    
    for (i=0;i<newSDE->N;i++)
    {
        int j;
        for (j=0;j<SSAL_MAX_NAME_SIZE-1;j++)
        {
            newSDE->names[i][j] = names[i][j];
        }
        newSDE->names[i][SSAL_MAX_NAME_SIZE-1] = '\0';
    }

    /*all done, return the new model*/
    return newModel;
}

/**
 * @brief write chemical reation network too file
 * @param stream the output stream
 * @param data a chemical reaction network
 *
 */
int SSAL_WriteChemicalReactionNetwork(FILE * stream ,SSAL_ChemicalReactionNetwork data)
{
    int i,j;

    if (stream == NULL)
    {
        return SSAL_IO_ERROR;
    }
    fprintf(stream,"{\"version\" : \"1\", \"level\" : \"3\",");
    fprintf(stream,"\"model\" : { \"id\" : \"crn\",\"name\" : \"crn\",\"substanceUnits\" : \"item\", \"timeUnits\" : \"seconds\",\"extentUnits\" : \"item\",");
    fprintf(stream,"\"compartments\" : [");
    fprintf(stream,"{\"id\" : \"space\",\"spatialDimensions\" : 0, \"size\" : 0}");
    fprintf(stream,"],");
    fprintf(stream,"\"species\" : [");
        fprintf(stream,"{\"id\" : \"%s\", \"compartment\" : \"space\",\"initialAmount\" : %d, \"hasOnlySubstanceUnits\" : true}",data.names[0],0);
    for (i=1;i<data.N;i++)
    {
        fprintf(stream,",{\"id\" : \"%s\", \"compartment\" : \"space\",\"initialAmount\" : %d, \"hasOnlySubstanceUnits\" : true}",data.names[i],0);
    }
    fprintf(stream,"],");
    fprintf(stream,"\"parameters\" : [");
        fprintf(stream,"{\"id\" : \"c%d\",\"value\" : %f}",0,data.c[0]);
    for (j=1;j<data.M;j++)
    {
        fprintf(stream,",{\"id\" : \"c%d\",\"value\" : %f}",j,data.c[j]);
    }
    fprintf(stream,"],");
    fprintf(stream,"\"reactions\" : [");
        fprintf(stream,"{\"id\" : \"R%d\",\"reversible\" : false, \"reactants\" : [",0);
        int first = 1;
        for (i=0;i<data.N;i++)
        {
            if (data.nu_minus[i] != 0)
            {
                if (first) first = 0;
                else fprintf(stream,",");
                fprintf(stream,"{\"species\" : \"%s\", \"stoichiometry\" : %.0f}",data.names[i],data.nu_minus[i]);
            }
        }
        fprintf(stream,"],\"products\" : [");
        first = 1;
        for (i=0;i<data.N;i++)
        {
            if (data.nu_plus[i] != 0)
            {
                if (first) first = 0;
                else fprintf(stream,",");

            fprintf(stream,"{\"species\" : \"%s\", \"stoichiometry\" : %.0f}",data.names[i],data.nu_plus[i]);
            }
        }
        fprintf(stream,"]}");
    for (j=1;j<data.M;j++)
    {
        fprintf(stream,",{\"id\" : \"R%d\",\"reversible\" : false, \"reactants\" : [",j);
        first = 1;
        for (i=0;i<data.N;i++)
        {
            if (data.nu_minus[j*data.N + i] != 0)
            {
                if (first) first = 0;
                else fprintf(stream,",");
            fprintf(stream,"{\"species\" : \"%s\", \"stoichiometry\" : %.0f}",data.names[i],data.nu_minus[j*data.N + i]);
            }
        }
        fprintf(stream,"],\"products\" : [");
        first = 1;
        for (i=0;i<data.N;i++)
        {
            if (data.nu_plus[j*data.N + i] != 0)
            {
                if (first) first = 0;
                else fprintf(stream,",");
            fprintf(stream,"{\"species\" : \"%s\", \"stoichiometry\" : %.0f}",data.names[i],data.nu_plus[j*data.N + i]);
            }
        }
        fprintf(stream,"]}");

    }
    fprintf(stream,"]");
    fprintf(stream,"}");
    fprintf(stream,"}\n");
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

    fprintf(stream,"\"T\"");
    for (j=0;j<sim->Nvar;j++)
    {
        fprintf(stream,",\"E[%s(T)]\"",sim->var[j]);
    }
    for (j=0;j<sim->Nvar;j++)
    {
        fprintf(stream,",\"V[%s(T)]\"",sim->var[j]);
    }
    fprintf(stream,"\n");
    /*print initial conditions*/
    fprintf(stream,"%f",0.0);
    for (j=0;j<sim->Nvar;j++)
    {
        fprintf(stream,",%f",sim->IC[sim->varInd[j]]); 
    }
    for (j=0;j<sim->Nvar;j++)
    {
       fprintf(stream,",%f",0.0); 
    }
    fprintf(stream,"\n");

    for (i=0;i<sim->NT;i++)
    {
        fprintf(stream,"%f",sim->T[i]);
        for (j=0;j<sim->Nvar;j++)
        {
            fprintf(stream,",%f",sim->E[j*sim->NT + i]);
        }
        for (j=0;j<sim->Nvar;j++)
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

    fprintf(stream,"\"T\",\"R\"");
    for (j=0;j<sim->Nvar;j++)
    {
        fprintf(stream,",\"%s(T)\"",sim->var[j]);
    }
    fprintf(stream,"\n");
    
    for (r=0;r<sim->NR;r++)
    {
        /*print initial conditions*/
        fprintf(stream,"%f,%d",0.0,r);
        for (j=0;j<sim->Nvar;j++)
        {
            fprintf(stream,",%f",sim->IC[sim->varInd[j]]); 
        }
        fprintf(stream,"\n");
        for (i = 0;i<sim->NT;i++)
        {
            fprintf(stream,"%f,%d",sim->T[i],r);
            for (j=0;j<sim->Nvar;j++)
            {
                fprintf(stream,",%f",sim->output[r*(sim->Nvar*sim->NT) + j*(sim->NT) + i]);
            }
            fprintf(stream,"\n");
        }
    }
}

/**
 * @brief imports chemical reaction network from an LSBML file.
 * @param filename the name of the lsbml file
 * @return A Chemical Reaction Network
 */
SSAL_Model SSAL_ImportLSBML(const char * filename )
{
    cJSON *root, *model, *parameters, *reactions, *species;
    FILE *fp;
    long int size;
    char * lsbml_json;
    int i,j;
    SSAL_ChemicalReactionNetwork *CRN;
    SSAL_Model newModel; 
    char * funcname;
    funcname = "SSAL_ImportLSBML";
    /*build wrapper struct*/
    CRN = (SSAL_ChemicalReactionNetwork*)malloc(sizeof(SSAL_ChemicalReactionNetwork));
    
    newModel.type = SSAL_CHEMICAL_REACTION_NETWORK;
    newModel.model = (void *)CRN;

    /*open the file*/
    if (!(fp = fopen(filename,"r")))
    {
        SSAL_HandleError(SSAL_IO_ERROR,funcname,__LINE__,1,0,NULL);
    }
    /*get file size in bytes*/
    fseek(fp,0,SEEK_END);
    size = ftell(fp);
    rewind(fp);
    /*allocate buffer*/
    if(!(lsbml_json = (char *)malloc(size*sizeof(char))))
    {
        SSAL_HandleError(SSAL_IO_ERROR,funcname,__LINE__,1,0,NULL);
    }
    /*read contents*/
    fread((void *)lsbml_json,sizeof(char),size,fp);

    /*parse the file*/
    root = cJSON_Parse(lsbml_json);
    model = cJSON_GetObjectItem(root,"model");
    parameters = cJSON_GetObjectItem(model,"parameters");
    reactions = cJSON_GetObjectItem(model,"reactions");
    species = cJSON_GetObjectItem(model,"species");

    CRN->M = cJSON_GetArraySize(reactions);
    CRN->N = cJSON_GetArraySize(species);  

    /*allocate memory for CRN components*/
    CRN->names = (char **)malloc((CRN->N)*sizeof(char*));
    for (i=0;i<CRN->N;i++)
    {
        CRN->names[i] = (char *)malloc(SSAL_MAX_NAME_SIZE*sizeof(char));
    }
    CRN->X0 = (SSAL_real_t*) malloc((CRN->N)*sizeof(SSAL_real_t));
    CRN->c = (SSAL_real_t*) malloc((CRN->M)*sizeof(SSAL_real_t));
    CRN->nu_minus = (SSAL_real_t*)malloc((CRN->N)*(CRN->M)*sizeof(SSAL_real_t));
    CRN->nu_plus = (SSAL_real_t*)malloc((CRN->N)*(CRN->M)*sizeof(SSAL_real_t));
    memset(CRN->nu_minus,0,(CRN->N)*(CRN->M)*sizeof(SSAL_real_t));
    memset(CRN->nu_plus,0,(CRN->N)*(CRN->M)*sizeof(SSAL_real_t));
    /*read species names*/
    for (i=0;i<CRN->N;i++)
    {
        strncpy(CRN->names[i],cJSON_GetObjectItem(cJSON_GetArrayItem(species,i),"id")->valuestring,SSAL_MAX_NAME_SIZE);
        CRN->X0[i] = (SSAL_real_t)(cJSON_GetObjectItem(cJSON_GetArrayItem(species,i),"initialAmount")->valuedouble);
    }
    /*read rate parameters*/
    /**@note we assume that the reaction rates are the only parameters, and that they are listed in the same order as the matching reaction*/
    for (j=0;j<CRN->M;j++)
    {
        CRN->c[j] = (SSAL_real_t)(cJSON_GetObjectItem(cJSON_GetArrayItem(parameters,j),"value")->valuedouble);
    }

    /* the reactions*/
    for (j=0;j<CRN->M;j++)
    {
       cJSON *react, *prod;
       int nReact, nProd,k;
       char temp[SSAL_MAX_NAME_SIZE];
       react = cJSON_GetObjectItem(cJSON_GetArrayItem(reactions,j),"reactants");
       prod = cJSON_GetObjectItem(cJSON_GetArrayItem(reactions,j),"products");

        nReact = cJSON_GetArraySize(react);
        nProd = cJSON_GetArraySize(prod);
        for (k=0;k<nReact;k++)
        {
            for (i=0;i<CRN->N;i++)
            {
                if (!strncmp(CRN->names[i],cJSON_GetObjectItem(cJSON_GetArrayItem(react,k),"species")->valuestring,SSAL_MAX_NAME_SIZE))
                {
                    CRN->nu_minus[j*(CRN->N)+i] = (SSAL_real_t)cJSON_GetObjectItem(cJSON_GetArrayItem(react,k),"stoichiometry")->valuedouble;
                }
            }
        }
        
        for (k=0;k<nProd;k++)
        {
            for (i=0;i<CRN->N;i++)
            {
                if (!strncmp(CRN->names[i],cJSON_GetObjectItem(cJSON_GetArrayItem(prod,k),"species")->valuestring,SSAL_MAX_NAME_SIZE))
                {
                    CRN->nu_plus[j*(CRN->N)+i] = (SSAL_real_t)cJSON_GetObjectItem(cJSON_GetArrayItem(prod,k),"stoichiometry")->valuedouble;
                }
            }
        }

    }
    cJSON_Delete(root);
    return newModel;
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
                                SSAL_real_t* T, SSAL_real_t * initCond)
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
    
    if((newRS->T = (SSAL_real_t *)malloc(NT*sizeof(SSAL_real_t)))==NULL)
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
    if((newRS->varInd = (int *)malloc(N*sizeof(int)))==NULL)
    {
        SSAL_HandleError(SSAL_IO_ERROR,funcname,__LINE__,1,0,NULL);
    }

    /*copy data */
    for (i=0;i<newRS->NT;i++)
    {
        newRS->T[i] = T[i];
    }
    for (i=0;i<newRS->Nvar;i++)
    {
        int j;
        if (obs != NULL)
        {
            strncpy(newRS->var[i],obs[i],SSAL_MAX_NAME_SIZE);
            newRS->var[i][SSAL_MAX_NAME_SIZE-1] = '\0';
        }
    } 

    /*setting initial condition sizes are dependent on the model*/
    switch (model->type)
    {
        case SSAL_CHEMICAL_REACTION_NETWORK:
        {
            SSAL_ChemicalReactionNetwork *CRN;
            CRN = (SSAL_ChemicalReactionNetwork *)model->model;
            SSAL_VARS2INDS(newRS,CRN,newRS->varInd)        
            /*initial conditions will be the initial chemical species copy numbers*/
            newRS->IC = (SSAL_real_t *)malloc((CRN->N)*sizeof(SSAL_real_t));
            for (i=0; i<(CRN->N); i++)
            {
                newRS->IC[i] = initCond[i];
            }

            newRS->output = (SSAL_real_t *)malloc((newRS->Nvar)*(newRS->NT)*(newRS->NR)*sizeof(SSAL_real_t));
        }
            break;
        case SSAL_STOCHASTIC_DIFFERENTIAL_EQUATION:
        {
            SSAL_StochasticDifferentialEquation *SDE;
            SDE = (SSAL_StochasticDifferentialEquation *)model->model;
            SSAL_VARS2INDS(newRS,SDE,newRS->varInd)        
            /*initial conditions will be the initial */
            newRS->IC = (SSAL_real_t *)malloc((SDE->N)*sizeof(SSAL_real_t));
            for (i=0; i<(SDE->N); i++)
            {
                newRS->IC[i] = initCond[i];
            }
            newRS->output = (SSAL_real_t *)malloc((newRS->Nvar)*(newRS->NT)*(newRS->NR)*sizeof(SSAL_real_t));
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
                                SSAL_real_t* T, SSAL_real_t * initCond)
{
    int i;
    char * obsArray;
    char * funcname;
    funcname = "SSAL_CreateExpectedValueSim";
    SSAL_Simulation newSim;
    SSAL_ExpectedValueSimulation *newEVS;
    newSim.type = SSAL_EXPECTEDVALUE;
    newSim.model = model;
    newEVS = (SSAL_ExpectedValueSimulation*)malloc(sizeof(SSAL_ExpectedValueSimulation));
    newSim.sim = (void*)newEVS;
    /*allocate memory*/
    newEVS->NR = NR;
    newEVS->NT = NT;
    newEVS->T = (SSAL_real_t *)malloc(NT*sizeof(SSAL_real_t));

    /**@todo extend to handle multiple scenarios */
    newEVS->NIC = 1; 
    newEVS->Nvar = N;
    obsArray = (char *)malloc(N*SSAL_MAX_NAME_SIZE*sizeof(char));
    newEVS->var = (char **)malloc(N*sizeof(char*));
    for (i=0;i<newEVS->Nvar;i++)
    {
        newEVS->var[i] = obsArray + i*SSAL_MAX_NAME_SIZE;
    }

    if((newEVS->varInd = (int *)malloc(N*sizeof(int)))==NULL)
    {
        SSAL_HandleError(SSAL_IO_ERROR,funcname,__LINE__,1,0,NULL);
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
        if (obs != NULL)
        {
            strncpy(newEVS->var[i],obs[i],SSAL_MAX_NAME_SIZE);
            newEVS->var[i][SSAL_MAX_NAME_SIZE-1] = '\0';
        }
    }


    /*setting initial condition sizes are dependent on the model*/
    switch (model->type)
    {
        case SSAL_CHEMICAL_REACTION_NETWORK:
        {
            SSAL_ChemicalReactionNetwork *CRN;
            CRN = (SSAL_ChemicalReactionNetwork *)model->model;
            SSAL_VARS2INDS(newEVS,CRN,newEVS->varInd) 
            /*initial conditions will be the initial chemical species copy numbers*/
            newEVS->IC = (SSAL_real_t *)malloc((CRN->N)*sizeof(SSAL_real_t));
            for (i=0; i<CRN->N; i++)
            {
                newEVS->IC[i] = initCond[i];
            }
            newEVS->E = (SSAL_real_t *)malloc((newEVS->Nvar)*(newEVS->NT)*sizeof(SSAL_real_t));
            newEVS->V = (SSAL_real_t *)malloc((newEVS->Nvar)*(newEVS->NT)*sizeof(SSAL_real_t));

        }
            break;
        case SSAL_STOCHASTIC_DIFFERENTIAL_EQUATION:
        {
            SSAL_StochasticDifferentialEquation *SDE;
            SDE = (SSAL_StochasticDifferentialEquation *)model->model;
            SSAL_VARS2INDS(newEVS,SDE,newEVS->varInd)        
            /*initial conditions will be the initial value*/
            newEVS->IC = (SSAL_real_t *)malloc((SDE->N)*sizeof(SSAL_real_t));
            for (i=0; i<(SDE->N); i++)
            {
                newEVS->IC[i] = initCond[i];
            }
            newEVS->E = (SSAL_real_t *)malloc((newEVS->Nvar)*(newEVS->NT)*sizeof(SSAL_real_t));
            newEVS->V = (SSAL_real_t *)malloc((newEVS->Nvar)*(newEVS->NT)*sizeof(SSAL_real_t));
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
        case SSAL_STOCHASTIC_DIFFERENTIAL_EQUATION:
            rc = SSAL_SimulateSDE(sim,alg,args);
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
    if (argv != NULL)
    {
        free(argv[0]);
        free(argv);
    }
    return rc;
}

/**
 * @brief Run a realisation simulation
 * @param sim a Realisation simulation structure
 * @param alg the selected algorithm type
 * @param args algorithm specific args
 */
int SSAL_SimulateSDE(SSAL_Simulation *sim, SSAL_AlgorithmType alg, 
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
            rc = SSAL_SimulateSDERealisations(sim->sim,sim->model->model,alg,argc,argv);
            break;
        case SSAL_EXPECTEDVALUE:
            rc = SSAL_SimulateSDEExpectedValue(sim->sim,sim->model->model,alg,argc,argv);
            break;
        default:
            
            break;
    }
    /*the first pointer is a pointer to the whole array*/
    if (argv != NULL)
    {
        free(argv[0]);
        free(argv);
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
    SSAL_real_t nu[model->N*model->M];
    SSAL_real_t * X_rj;

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
#ifdef __FLOAT64__
                degils(model->M,model->N,sim->NT,sim->T,sim->IC,model->nu_minus,
                    nu,model->c,sim->Nvar,sim->varInd,X_rj+j*(sim->NT*sim->Nvar));
#else
                segils(model->M,model->N,sim->NT,sim->T,sim->IC,model->nu_minus,
                    nu,model->c,sim->Nvar,sim->varInd,X_rj+j*(sim->NT*sim->Nvar));
#endif
            }
        }
            break;
        case SSAL_ASSA_EULER_MARUYAMA_SEQUENTIAL:
        case SSAL_ASSA_TAU_LEAP_SEQUENTIAL:
        {
            SSAL_real_t tau;
            tau = (SSAL_real_t)atof(SSAL_GetArg("--tau",argc,argv));

            for (j=0;j<sim->NR;j++)
            {
#ifdef __FLOAT64__
                datauls(model->M,model->N,sim->NT,sim->T,sim->IC,model->nu_minus,
                    nu,model->c,sim->Nvar,sim->varInd,tau,X_rj+j*(sim->NT*sim->Nvar));
#else
                satauls(model->M,model->N,sim->NT,sim->T,sim->IC,model->nu_minus,
                    nu,model->c,sim->Nvar,sim->varInd,tau,X_rj+j*(sim->NT*sim->Nvar));
#endif
            }
        }
            break;
        default:
            SSAL_HandleError(SSAL_UNSUPPORTED_ALGORITHM_ERROR,"SSAL_SimulateCRN",__LINE__,1,0,NULL);
            break;
    }

}

/**
 * @brief run a realisation simulation on a SDE model
 * @param sim a Realisation Simulation structure
 * @param model an SDE model
 * @param alg the selected algorithm type
 * @param argc the number of args
 * @param argv input ags
 */
int SSAL_SimulateSDERealisations(SSAL_RealisationSimulation *sim, 
            SSAL_StochasticDifferentialEquation *model, SSAL_AlgorithmType alg, int argc, char ** argv)
{
    int j,i;
    SSAL_real_t * X_rj;

    X_rj = sim->output;
    /*algorithm selector*/ 
    switch(alg)
    {
        case SSAL_ASSA_EULER_MARUYAMA_SEQUENTIAL:
        {
            SSAL_real_t h;
            h = (SSAL_real_t)atof(SSAL_GetArg("--h",argc,argv));
            for (j=0;j<sim->NR;j++)
            {
#ifdef __FLOAT64__
                daems(model->M,model->N,sim->NT,sim->T,model->p,sim->IC,model->mu,model->sigma
                    ,sim->Nvar,sim->varInd,h,X_rj+j*(sim->NT*sim->Nvar));
#else
                saems(model->M,model->N,sim->NT,sim->T,model->p,sim->IC,model->mu,model->sigma
                    ,sim->Nvar,sim->varInd,h,X_rj+j*(sim->NT*sim->Nvar));
#endif
            }
        }
            break;
        case SSAL_ASSA_EULER_MARUYAMA_CORRELATED_SEQUENTIAL:
        {
            SSAL_real_t h_f;
            int32_t M;
            h_f = (SSAL_real_t)atof(SSAL_GetArg("--h",argc,argv));
            M = (int32_t)atoi(SSAL_GetArg("--M",argc,argv));
            if ((sim->NR)%2 != 0)
            {
                    SSAL_HandleError(SSAL_UNSUPPORTED_ALGORITHM_ERROR,"SSAL_SimulateCRN",__LINE__,1,0,NULL);
                    return 1;
            }
            else
            {
                for (j=0;j<sim->NR;j+=2)
                {

#ifdef __FLOAT64__
                    dacems(model->M,model->N,sim->NT,sim->T,model->p,sim->IC,model->mu,model->sigma,sim->Nvar,sim->varInd,h_f,M,X_rj + j*(sim->NT*sim->Nvar),X_rj + (j+1)*(sim->NT*sim->Nvar));
#else
                    SSAL_HandleError(SSAL_UNSUPPORTED_ALGORITHM_ERROR,"SSAL_SimulateCRN",__LINE__,1,0,NULL);
                    return 1;
#endif
                }
            }
        }
            break;
        default:
            SSAL_HandleError(SSAL_UNSUPPORTED_ALGORITHM_ERROR,"SSAL_SimulateCRN",__LINE__,1,0,NULL);
            return 1;
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
    SSAL_real_t nu[model->N*model->M];
    SSAL_real_t * X_r;
    SSAL_real_t * E_X;
    SSAL_real_t * V_X;

    for (i=0;i<model->N*model->M;i++)
    {
        nu[i] = model->nu_plus[i] - model->nu_minus[i] ;
    }

    X_r = (SSAL_real_t *)malloc(sim->NT*sim->Nvar*sizeof(SSAL_real_t));
    E_X = (SSAL_real_t *)malloc(sim->NT*sim->Nvar*sizeof(SSAL_real_t));
    V_X = (SSAL_real_t *)malloc(sim->NT*sim->Nvar*sizeof(SSAL_real_t));
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        E_X[i] = 0;
    }
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        V_X[i] = 0;
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
#ifdef __FLOAT64__
                degils(model->M,model->N,sim->NT,sim->T,sim->IC,
                    model->nu_minus,nu,model->c,sim->Nvar,sim->varInd,X_r);
#else
                segils(model->M,model->N,sim->NT,sim->T,sim->IC,
                    model->nu_minus,nu,model->c,sim->Nvar,sim->varInd,X_r);
#endif
                
                for (i=0;i<sim->Nvar*sim->NT;i++)
                {
                    E_X[i] += X_r[i];
                }
                for (i=0;i<sim->Nvar*sim->NT;i++)
                {
                    V_X[i] += X_r[i]*X_r[i];
                }
            }

            for (i=0;i<sim->Nvar*sim->NT;i++)
            {
                E_X[i] /= (SSAL_real_t)(sim->NR);
            }
            for (i=0;i<sim->Nvar*sim->NT;i++)
            {
                V_X[i] /= (SSAL_real_t)(sim->NR);
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
                V_X[i] /= (SSAL_real_t)(sim->NR);
            }

        }
            break;
        case SSAL_ASSA_EULER_MARUYAMA_SEQUENTIAL:
        case SSAL_ASSA_TAU_LEAP_SEQUENTIAL:
        {
            SSAL_real_t tau;
            tau = (SSAL_real_t)atof(SSAL_GetArg("--tau",argc,argv));
            
            for (j=0;j<sim->NR;j++)
            {
#ifdef __FLOAT64__
                datauls(model->M,model->N,sim->NT,sim->T,sim->IC,
                    model->nu_minus,nu,model->c,sim->Nvar,sim->varInd,tau,X_r);
#else
                satauls(model->M,model->N,sim->NT,sim->T,sim->IC,
                    model->nu_minus,nu,model->c,sim->Nvar,sim->varInd,tau,X_r);
#endif
                for (i=0;i<sim->Nvar*sim->NT;i++)
                {
                    E_X[i] += (double)X_r[i];
                }
                for (i=0;i<sim->Nvar*sim->NT;i++)
                {
                    V_X[i] +=(double)( X_r[i]*X_r[i]);
                }
            }

            for (i=0;i<sim->Nvar*sim->NT;i++)
            {
                E_X[i] /= (SSAL_real_t)(sim->NR);
            }
            for (i=0;i<sim->Nvar*sim->NT;i++)
            {
                V_X[i] /= (SSAL_real_t)(sim->NR);
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
                V_X[i] /= (SSAL_real_t)(sim->NR);
            }


        }
            break;
        case SSAL_ASSA_MULTI_LEVEL_SEQUENTIAL:
        {
            SSAL_real_t tau0, eps;
            int M, L;
            tau0 = (SSAL_real_t)atof(SSAL_GetArg("--tau0",argc,argv));
            M = (int)atoi(SSAL_GetArg("--M",argc,argv));
            L = (int)atoi(SSAL_GetArg("--L",argc,argv));
            eps = (SSAL_real_t)atof(SSAL_GetArg("--eps",argc,argv));
#ifdef __FLOAT64__
            damlmcbs(model->M,model->N,sim->NT,sim->T,sim->IC,model->nu_minus,
                nu,model->c,tau0,M,L,eps,sim->Nvar,sim->varInd,NULL,E_X,V_X);
#else
            samlmcbs(model->M,model->N,sim->NT,sim->T,sim->IC,model->nu_minus,
                nu,model->c,tau0,M,L,eps,sim->Nvar,sim->varInd,NULL,E_X,V_X);
#endif
        }
            break;
        default:
            SSAL_HandleError(SSAL_UNSUPPORTED_ALGORITHM_ERROR,"SSAL_SimulateCRN",__LINE__,1,0,NULL);
            break;
    }

    
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        sim->E[i] =  E_X[i];
        sim->V[i] =  V_X[i];
    }

    free(X_r);
    free(E_X);
    free(V_X);
}
 
 /**
 * @brief Compute Expected value on a SDE model
 * @param sim an Expected Value Simulation structure
 * @param model a SDE model
 * @param alg the selected algorithm type
 * @param argc the number of args
 * @param argv input ags
 */
int SSAL_SimulateSDEExpectedValue(SSAL_ExpectedValueSimulation *sim, 
            SSAL_StochasticDifferentialEquation *model, SSAL_AlgorithmType alg, int argc, char ** argv)
{
    int j,i;
    SSAL_real_t * X_r;
    SSAL_real_t * E_X;
    SSAL_real_t * V_X;


    X_r = (SSAL_real_t *)malloc(sim->NT*sim->Nvar*sizeof(SSAL_real_t));
    E_X = (SSAL_real_t *)malloc(sim->NT*sim->Nvar*sizeof(SSAL_real_t));
    V_X = (SSAL_real_t *)malloc(sim->NT*sim->Nvar*sizeof(SSAL_real_t));
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        E_X[i] = 0;
    }
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        V_X[i] = 0;
    }
    /*algorithm selector*/ 
    switch(alg)
    {
        case SSAL_ASSA_EULER_MARUYAMA_SEQUENTIAL:
        {
            SSAL_real_t h;
            h = (SSAL_real_t)atof(SSAL_GetArg("--h",argc,argv));
            
            for (j=0;j<sim->NR;j++)
            {
#ifdef __FLOAT64__
                daems(model->M,model->N,sim->NT,sim->T,model->p,sim->IC,
                    model->mu,model->sigma,sim->Nvar,sim->varInd,h,X_r);
#else
                saems(model->M,model->N,sim->NT,sim->T,model->p,sim->IC,
                    model->mu,model->sigma,sim->Nvar,sim->varInd,h,X_r);
#endif
                for (i=0;i<sim->Nvar*sim->NT;i++)
                {
                    E_X[i] += (double)X_r[i];
                }
                for (i=0;i<sim->Nvar*sim->NT;i++)
                {
                    V_X[i] +=(double)( X_r[i]*X_r[i]);
                }
            }

            for (i=0;i<sim->Nvar*sim->NT;i++)
            {
                E_X[i] /= (SSAL_real_t)(sim->NR);
            }
            for (i=0;i<sim->Nvar*sim->NT;i++)
            {
                V_X[i] /= (SSAL_real_t)(sim->NR);
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
                V_X[i] /= (SSAL_real_t)(sim->NR);
            }


        }
            break;
        default:
            SSAL_HandleError(SSAL_UNSUPPORTED_ALGORITHM_ERROR,"SSAL_SimulateCRN",__LINE__,1,0,NULL);
            break;
    }

    
    for (i=0;i<sim->Nvar*sim->NT;i++)
    {
        sim->E[i] =  E_X[i];
        sim->V[i] =  V_X[i];
    }

    free(X_r);
    free(E_X);
    free(V_X);
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
    if (args[0] == '\0')
    {
        *argc = 0;
        return NULL;
    }
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

