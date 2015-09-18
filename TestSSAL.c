#include "SSAL.h"

int main(int argc , char ** argv)
{
    int i;
    int n,m;
    int N;
    int NT = 5;
    float dt;
    float T[3];
    n = 3;
    m = 5;
    N = 100;
    /*this is just the Dimerisation problem*/
    float nu_minus[5*3] = {0,0,0,1,0,0,0,2,0,1,0,0,0,1,0};
    float nu_plus[5*3] = {1,0,0,1,1,0,0,0,1,0,0,0,0,0,0};
    float c[5] = {25,1000,0.001,0.1,1};
    char *names[3] = {"M","P","D"};
    float X0[3] = {0,0,0};
    SSAL_Model CRN;
    SSAL_Simulation sim;

    /*build sample times*/
    T[NT-1] = 1.0;
    dt = 1.0/((float)(NT));
    for (i=0; i < NT-1; i++){
        T[i] = dt*((float)(i+1));
    }

    srand(1337);
    /*init SSAL */
    SSAL_Initialise(argc,argv);
    /* build chemical reaction network*/
    CRN = SSAL_CreateChemicalReactionNetwork(names,m,n,nu_minus,nu_plus,c);
    /* build realisation simulation */
    //sim = SSAL_CreateRealisationsSim(&CRN,n,names,N,NT,T,X0);
    sim = SSAL_CreateExpectedValueSim(&CRN,n,names,N,NT,T,X0);
    /*write the CRN*/
    //SSAL_WriteChemicalReactionNetwork(stdout,*((SSAL_ChemicalReactionNetwork *)(CRN.model)));
    /*simulate realisations*/
    SSAL_Simulate(&sim,SSAL_ESSA_GILLESPIE_SEQUENTIAL,"");
    /*write the data*/
    //SSAL_WriteRealisationsSim(stdout,sim.sim);
    return 0;
}
