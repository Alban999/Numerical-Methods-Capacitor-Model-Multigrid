#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "GS.h"

int GS(int **ia, int **ja, double **a, double **b, double **r, int preOrPost, double **u, int n)
/*
 
 But
 ===
 Effectuer la méthode de Gauss Seidel
 
 Arguments
 =========
 ia (output) - pointeur vers le tableau 'ia' de la matrice A
 ja (output) - pointeur vers le tableau 'ja' de la matrice A
 a  (output) - pointeur vers le tableau 'a' de la matrice A
 b  (output) - pointeur vers le tableau 'b'
 r  (output) - pointeur vers le tableau 'r' (residu)
 preOrPost  (input) - choix de pre- (1) or post-smoothing
 u  (output) - pointeur vers le tableau 'u' (solution)
 n  (input) - pointeur vers le nombre d'inconnus dans le système
 
 */
{
    
    /* DECLARATION DES VARIABLES */
    
    int i = 0, j = 0, k = 0;
    double termDiag = 0, tot = 0, tot1 = 0;
    
    //preOrPost = 1 -> pré smoothing | preOrPost = 0 -> post smoothing
    if(preOrPost){
        for(i = 0; i<n; i++){
            (*u)[i] = (*b)[i];
            
            for(j = (*ia)[i]; j < (*ia)[i+1]; j++){
                
                if((*ja)[j]==i){
                    
                    termDiag = (*a)[j];
                }
                else{
                    
                    (*u)[i]  -= (*a)[j]*(*u)[(*ja)[j]];
                }
            }
            
            (*u)[i] /= termDiag;
            
        }
    }
    else{
        for(i = n-1; i>-1; i--){
            (*u)[i]  = (*b)[i];
            
            for(j = (*ia)[i]; j < (*ia)[i+1]; j++){
                if((*ja)[j]==i){
                    termDiag = (*a)[j];
                }
                else{
                    (*u)[i]  -= (*a)[j]*(*u)[(*ja)[j]];
                }
            }
            
            (*u)[i] /= termDiag;
        }
    }
    
    /* CALCUL DU RESIDU */
    
    for(i = 0; i< n; i++){
        (*r)[i] = (*b)[i];
        
        
        for(k = (*ia)[i]; k<(*ia)[i+1]; k++)
        {
            
            (*r)[i] -= (*a)[k]*(*u)[(*ja)[k]];
            
            
        }
        
        tot += (*r)[i]*(*r)[i];
        tot1 += (*b)[i]*(*b)[i];
    }
    
    if(!preOrPost){
        //Affichage de ||r||/||b||
        
        //printf("\nresidu %e\n", sqrt(tot)/sqrt(tot1));
    }
    
    
    return 0;
}

