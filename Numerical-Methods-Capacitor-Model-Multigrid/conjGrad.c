#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "residu.h"
#include "multigrille.h"
#include "conjGrad.h"

int conjGrad(int ***ia, int ***ja, double ***a, double ***b, double **u, double ***uTab, int **nTab, double h, double bornes[4], double Lx, double Ly, int *niveau, double hcMaxGrossier, int nbreIter, double **rNormMG, int nbreIterMax)
/*
 
 But
 ===
 Résoudre le système Au = b à l'aide du gradient conjugué
 
 Arguments
 =========
 iaTab (input)  - tableau de pointeurs 'ia'
 aTab  (input)  - tableau de pointeurs 'a'
 jaTab (input)  - tableau de pointeurs 'ja'
 bTab  (input)  - tableau de pointeurs 'b'
 uTab  (output) - tableau de pointeurs correspondant chacun au solution de chaque niveau
 nTab  (input)  - tableau contenant le nombre d'inconnues des différents niveaux
 h     (input)  - le pas de discrétisation
 bornes  (input) - les bornes entre lesquelles est défini le condensateur
 Lx      (input) - longeur du condensateur dans la direction 1x
 Ly      (input) - longeur du condensateur dans la direction 1y
 *niveau (input) - nombre correspondant au niveau auquel on se trouve
 hcMaxGrossier (input) - le pas le plus grossier (pas auquel on résout à l'aide de UMFPACK)
 nbreItMG (input) - numéro de l'itération
 **rNorm  (input) - tableau contenant les normes des résidus pour chaque itération
 r      (output) - pointeur vers le tableau 'r' (residu)
 rc     (output) - pointeur vers le tableau 'rc' (residu restreint)
 mxc    (input)  - nombre de points dans la direction 1x sur la grille
 myc    (input)  - nombre de points dans la direction 1y sur la grille
 bornes (input)  - les bornes entre lesquelles est défini le condensateur
 h      (input)  - le pas de discrétisation de la grille fine
 hc     (input)  - le pas de discrétisation de la grille grossière
 nx     (input)  - le nombre d'éléments à l'intérieur du condensateur dans la directions 1x
 ny     (input)  - le nombre d'éléments à l'intérieur du condensateur dans la directions 1x
 
 */
{
    /* DECLARATION DE VARIABLES */
    
    double beta = 0, rsSold = 0, rsNew = 0, alpha = 0, numAlpha = 0, denAlpha = 0, tot = 0, tot1 = 0;
    int i = 0, n = (*nTab)[0];
    double *Ad, *r, *d;
    
    
    //Inisialisation
    Ad = malloc(n*sizeof(double));
    r = malloc(n*sizeof(double));
    d = malloc(n*sizeof(double));
    
    //Initialisation
    for(i=0;i<n;i++){
        (*u)[i] = 0;
    }
    
    if(Ad == NULL || r == NULL || d == NULL || b == NULL){
        printf("ERREUR : problème d'allocation");
        return 1;
    }
    
    //Calcul de la norme de b
    for(i=0;i<n;i++){
        tot1 += (*b)[0][i]*(*b)[0][i];
    }
    
    resVect(&(*u), &(*a)[0], &(*ia)[0], &(*ja)[0], &(*b)[0], n, &r);
    
    for(i=0; i<n; i++){
        (*b)[0][i] = r[i];
    }
    
    (*niveau) = -1;
    multigrille(uTab, nTab, h, a, ia, ja, b, bornes, Lx, Ly, niveau, hcMaxGrossier, nbreIter, rNormMG);
    
    
    for(i=0; i<n; i++){
        d[i] = (*uTab)[0][i];
        Ad[i] = 0;
        rsSold += r[i]*(*uTab)[0][i];
    }
    
    
    
    
    //Commencement des itérations
    for(int iter = 1; iter<nbreIterMax-1; iter++){
        
        //Calcul de Ad
        for(int t=0; t<n; t++){
            Ad[t] = 0;
            for(int j = (*ia)[0][t]; j < (*ia)[0][t+1]; j++){
                Ad[t] += (*a)[0][j]*d[(*ja)[0][j]];
            }
        }
        
        //Calcul de alpha
        denAlpha = 0;
        for(int t=0; t<n; t++){
            denAlpha += d[t]*Ad[t];
        }
        
        alpha = rsSold/denAlpha;
        
        for(int t=0; t<n; t++){
            (*u)[t] += alpha*d[t];
            r[t] -= alpha*Ad[t];
            
        }
        
        
        for(i=0; i<n; i++){
            (*b)[0][i] = r[i];
        }
        
        (*niveau) = -1;
        multigrille(uTab, nTab, h, a, ia, ja, b, bornes, Lx, Ly, niveau, hcMaxGrossier, nbreIter, rNormMG);
        
        rsNew = 0;
        for(int t=0; t<n; t++){
            rsNew+= r[t]*(*uTab)[0][t];
        }
        
        //Calcul de beta
        beta = rsNew/rsSold;
        
        
        for(int t=0; t<n; t++){
            d[t] = (*uTab)[0][t] + beta*d[t];
        }
        
        tot = 0;
        for(int t=0; t<n; t++){
            tot += r[t]*r[t];
        }
        
        //Affichage du résidu
        printf("\nresidu %e\n", sqrt(tot)/sqrt(tot1));
        rsSold = rsNew;
        
        
        
    }
    
    
    return 0;
}

