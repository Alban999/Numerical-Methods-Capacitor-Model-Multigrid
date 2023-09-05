
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "residu.h"

int resVect(double **u, double **a, int **ia, int **ja, double **b, int n, double **r)
/*
 
 But
 ===
 Calculer le vecteur résidu du système linéaire Ax = b
 
 Arguments
 =========
 u  (input)  - le vecteur solution
 a  (input)  - pointeur vers le tableau 'a' de la m
 ia (input)  - pointeur vers le tableau 'ia' de la matrice A
 ja (input)  - pointeur vers le tableau 'ja' de la matrice A
 b  (input)  - pointeur vers le tableau 'b'
 n  (input)  - le nombre d'inconnus dans le système
 r  (output) - le vecteur résidu

 */
{
    
    for(int i = 0; i< n; i++){
        (*r)[i] = (*b)[i];
        
        
        for(int k = (*ia)[i]; k<(*ia)[i+1]; k++)
        {
            
            (*r)[i] -= (*a)[k]*(*u)[(*ja)[k]];
            
            
        }

    }
    
    return 0;
}

