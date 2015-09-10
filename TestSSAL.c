#include "SSAL.h"

int main()
{
    int n,m;
    n = 4;
    m = 5;
    float nu_minus[5*4] = {1,0,0,0,0,1,0,0,0,0,2,0,0,1,0,0,0,0,1,0};
    float nu_plus[5*4] = {1,1,0,0,0,1,1,0,0,0,0,1,0,0,0,0,0,0,0,0};
    float c[5] = {25,1000,0.001,0.1,1};
    char *names[4] = {"G","M","P","D"};
    SSAL_ChemicalReactionNetwork CRN;

    CRN = SSAL_CreateChemicalReactionNetwork(names,m,n,nu_minus,nu_plus,c);
    SSAL_WriteChemicalReactionNetwork(stdout,CRN);
}
