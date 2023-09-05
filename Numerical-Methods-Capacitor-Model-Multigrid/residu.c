#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "residu.h"

double residu(double **x, double **a, int **ia, int **ja, double **b, int n)
/*
 
 But
 ===
 Calculer le residu du système linéaire Ax = b grâce à la formule suivante:
 
            ||Ax-b||_2
 residu  =  --------
             ||b||_2
 
 Arguments
 =========
 x  (input) - le vecteur solution
 a  (input) - pointeur vers le tableau 'a' de la m
 ia (input) - pointeur vers le tableau 'ia' de la matrice A
 ja (input) - pointeur vers le tableau 'ja' de la matrice A
 b  (input) - pointeur vers le tableau 'b'
 n  (input) - le nombre d'inconnus dans le système
 
 */
{
    
    /* DECLARATION DES VARIABLES */
    
    int i = 0, j = 0;
    double total = 0, residu = 0, totb = 0;
    double ligne = 0;
    
    /* CALCUL DU RESIDU */
        
    for(i = 0; i<n; i++)
    {
        ligne = -(*b)[i];
        for(j = (*ia)[i]; j<(*ia)[i+1]; j++)
        {
            ligne += ((*a)[j]*(*x)[(*ja)[j]]);
        }
        total += ligne*ligne;
        totb += (*b)[i]*(*b)[i];
    }
    
    //residu = ||Ax-b||_2  / ||b||_2
    residu = sqrt(total)/sqrt(totb);

    return residu;
}

