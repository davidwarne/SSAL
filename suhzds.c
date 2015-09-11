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
int suhzds(int m,int n,float * restrict R, float * c, float * restrict X, float * restrict A)
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

