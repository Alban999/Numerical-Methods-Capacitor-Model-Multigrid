#include <stdio.h>
#include <stdlib.h>
#include "agmgSolveur.h"

void dagmg_(int*,double*,int*,int*,double*,double*,int*,int*,int*,int*,double*);

int solve_AGMG(int mx, int my, int n, double *a, int *ja, int *ia, double *b, double *x, int zero, int one, int iter, double tol, int iprint)
/*
 
 But
 ===
 Résoudre le système linéaire
 
 Ax = b
 
 en utilisant les fonctions de AGMG. La matrice A est
 supposée être dans le format CRS défini par les trois
 tableaux : 'ia', 'ja' et 'a'. Attention, ici le solveur
 est en FORTRAN. Il faut donc rajouter +1 aux tableaux
 'ia' et 'ja'.
 
 Arguments
 =========
 mx (input) - nombre de points dans la direction 1x sur la grille
 my (input) - nombre de points dans la direction 1y sur la grille
 n  (input) - le nombre d'inconnus dans le système
 a  (output) - tableau 'a' de la matrice A
 ja (output) - le tableau 'ja' de la matrice A
 ia (output) - le tableau 'ia' de la matrice A
 b  (output) - le tableau 'b'
 x  (output) - le vecteur solution
 zero (input) - 0
 one (input) - 1
 iter (input) - nombre maximal d'itérations
 tol (input) - tolérance sur la norme relative du résidu
 iprint (intput) - nombre pour les messages de sorties (6 = sortie standard)

 */
{
    /* DECLARATION DES VARIABLES */
    
    
    int i = 0;
    int *iaCopie, *jaCopie;
    double *bCopie;
    
    
    //Copi de ia et ja
    iaCopie  = malloc((n + 1) * sizeof(int));
    jaCopie  = malloc(ia[n] * sizeof(int));
    bCopie = malloc(n * sizeof(double));
    
    if ( iaCopie == NULL || jaCopie == NULL || bCopie == NULL) {
        printf("\n ERREUR : pas de mémoire pour les tableaux (AGMG)\n\n");
        return 1;
    }
    
    for(i=0;i<n+1;i++){
        iaCopie[i] = ia[i]+1;
    }
    
    for(i=0;i<ia[n];i++){
        jaCopie[i] = ja[i]+1;
    }
    
    for(i=0;i<n;i++){
        bCopie[i] = b[i];
    }
    
    dagmg_(&n,a,jaCopie,iaCopie,bCopie,x,&zero,&iprint,&one,&iter,&tol);
    
    free(jaCopie); free(iaCopie); free(bCopie);
                        
    return 0;
}
