#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "conjugateGrad.h"
#include "residu.h"
#include "multigrille.h"
#include "GC.h"

int CG(int **ia, int **ja, double **a, double **b, double **u, int n){
    /* DECLARATION DE VARIABLES */
    double beta = 0, rsSold = 0, rsNew = 0, alpha = 0, numAlpha = 0, denAlpha = 0, tot = 0, tot1 = 0;
    int i = 0;
    double *Ad, *r, *d;
    Ad = malloc(n*sizeof(double));
    r = malloc(n*sizeof(double));
    d = malloc(n*sizeof(double));
    
    if(Ad == NULL || r == NULL || d == NULL ){
        printf("ERREUR : problème d'allocation");
        return 1;
    }
    
    printf("\n\n\t- Displaying the matrix in CSR format -\n\n");
    
    i = 0;
    for(int line = 0; line < n; line++){
        for(int col = 0; col < n; col++){
            if (col == (*ja)[i]){
                printf("%f ", (*a)[i]);
                i++;
            } else {
                printf("%d ",0);
            }
        }
        printf("\n");
    }
    printf("------------------------------------\n");
    //Initialisation
    for(i=0;i<n;i++){
        (*u)[i] = 0;
    }
    
    if (solve_umfpack(n, (*ia), (*ja), (*a), (*b), (*u)))
        return 1;
    
    
    //Initialisation
    for(i=0;i<n;i++){
        (*u)[i] = 0;
    }
    
    resVect(&(*u), &(*a), &(*ia), &(*ja), &(*b), n, &r);
    
    for(i=0; i<n; i++){
        d[i] = r[i];
        Ad[i] = 0;
        rsSold += r[i]*r[i];
    }
    
    
    
    
    //Commencement des itérations
    for(int iter = 1; iter<n+1; iter++){
        //Calcul de Ad
        for(int t=0; t<n; t++){
            Ad[t] = 0;
            for(int j = (*ia)[t]; j < (*ia)[t+1]; j++){
                Ad[t] += (*a)[j]*d[(*ja)[j]];
            }
        }
        
        
        
        denAlpha = 0;
        for(int t=0; t<n; t++){
            denAlpha += d[t]*Ad[t];
        }
        
        alpha = rsSold/denAlpha;

        for(int t=0; t<n; t++){
            (*u)[t] = (*u)[t] + alpha*d[t];
            r[t] = r[t] - alpha*Ad[t];
            
        }
        
        printf("%.30e Y\n", r[1]);
        
        
        printf("\nres ICI\n");
        
        printf("%e %e \n", r[0], r[1]);
        
        rsNew = 0;
        for(int t=0; t<n; t++){
            rsNew+= r[t]*r[t];
        }
        
        
        
        
        beta = rsNew/rsSold;
        
        printf("\nrsSold %e\n", rsSold);
        printf("\nrsNew %e\n", rsNew);
        
        
        
        
        
        
        
        
        //printf("\Ad %f %f\n", Ad[0], Ad[1]);
        
        
        
        printf("\nbeta %e\n", beta);
        
        for(int t=0; t<n; t++){
            d[t] = r[t] + beta*d[t];
        }
        
        printf("\ndsol %e %e\n", d[0], d[1]);
        
        
        
        for(int t=0; t<n; t++){
            tot += r[t]*r[t];
            tot1 += (*b)[t]*(*b)[t];
        }
        
        printf("\nresidu %e\n", sqrt(tot)/sqrt(tot1));
        
        rsSold = rsNew;
        
        
        
    }
    printf("\nsSoL %f %f %f\n", (*u)[0], (*u)[1], (*u)[2]);
        
    
    return 0;
}
