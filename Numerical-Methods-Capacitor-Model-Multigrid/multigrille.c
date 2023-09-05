#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "multigrille.h"
#include "rho.h"
#include "prob.h"
#include "restriction.h"
#include "prolongation.h"
#include "GS.h"
#include "nombreElem.h"
#include "umfpk.h"
#include "tableauPotentiel.h"
#include "residu.h"
#include "resVect.h"

int multigrille(double ***uTab, int **nTab, double h, double ***aTab, int ***iaTab, int ***jaTab, double ***bTab, double bornes[4], double Lx, double Ly, int *niveau, double hcMaxGrossier, int nbreItMG, double **rNorm)
/*
 
 But
 ===
 Renvoyer la solution du système Au = b
 
 Arguments
 =========
 uTab (output) - tableau de pointeurs correspondant chacun au solution de chaque niveau
 nTab  (input) - tableau contenant le nombre d'inconnues des différents niveaux
 h       (input) - le pas de discrétisation
 aTab (input) - tableau de pointeurs 'a'
 iaTab (input) - tableau de pointeurs 'a'
 jaTab (input) - tableau de pointeurs 'a'
 bTab (input) - tableau de pointeurs 'a'
 bornes  (input) - les bornes entre lesquelles est défini le condensateur
 Lx      (input) - longeur du condensateur dans la direction 1x
 Ly      (input) - longeur du condensateur dans la direction 1y
 niveau (input) - nombre correspondant au niveau auquel on se trouve
 hcMaxGrossier (input) - le pas le plus grossier (pas auquel on résout à l'aide de UMFPACK)
 nbreItMG (input) - numéro de l'itération
 rNorm  (input) - tableau contenant les normes des résidus pour chaque itération
 
 */
{
    
    /* DECLARATION DES VARIABLES */
    
    //Dimensions du condensateur
    double hc = 2*h, tot1=0;
    int my = Ly/h+1, mx = Lx/h+1, mxc = Lx/hc+1, myc = Ly/hc+1;
    int niveauRel = 0;

    //Potentiel sur le bord bas du condensateur
    double uB = -0.5;

    //Autres variables
    int i;
    
    //Tableau à afficher graphiquement
    double *tableau;
    
    //Afficher les graphes de résidus : oui ou non
    int afficherRes = 0;
    
    //Niveau en plus
    (*niveau)++;

    //Stocker le niveau dans une varibale locale propre à un certain niveau
    niveauRel = *niveau;

    //Déclaration du vecteur résidu
    double *r, *rBis;
    r = malloc((*nTab)[niveauRel] * sizeof(double));
    rBis = malloc((*nTab)[niveauRel] * sizeof(double));
    
    if(r==NULL || rBis==NULL){
        printf("\nErreur : problème d'allocation\n");
        return 1;
    }
    
    /* DEBUT DE L'ALGORITHME */

    //Afficher résidu
    if(niveauRel==0){
        for(i=0; i<(*nTab)[0]; i++){
            tot1 += (*bTab)[0][i]*(*bTab)[0][i];
        }
        for(i=0; i<(*nTab)[0]; i++){
            r[i] = (*bTab)[0][i];
        }
        for(i=0; i<(*nTab)[0]; i++){
            rBis[i]=r[i]/sqrt(tot1);
        }
        if((nbreItMG==0 || nbreItMG==1) && afficherRes){
            printf("\n\nResidu avant pré-smoothing \n\n");
            tableau = malloc(mx*my*sizeof(double));
            creaTableau(rBis, tableau, mx, my, bornes, h, 0.0);
            if (plot(rBis, mx, my, h, tableau, "dataResidu.txt", "fonctionResidu.txt"))
                return 1;
        }
    }
    
    //Pré-smoothing
    for(i=0;i<4;i++){
        GS(&(*iaTab)[niveauRel], &(*jaTab)[niveauRel], &(*aTab)[niveauRel], &(*bTab)[niveauRel], &r, 1, &(*uTab)[niveauRel], (*nTab)[niveauRel]);
    }

    
    //Afficher résidu
    if(niveauRel==0){
        for(i=0; i<(*nTab)[0]; i++){
            rBis[i]=r[i]/sqrt(tot1);
        }
        if((nbreItMG==0 || nbreItMG==1) && afficherRes){
            printf("\n\nResidu après pré-smoothing \n\n");
            tableau = malloc(mx*my*sizeof(double));
            creaTableau(rBis, tableau, mx, my, bornes, h, 0.0);
            if (plot(rBis, mx, my, h, tableau, "dataResidu.txt", "fonctionResidu.txt"))
                return 1;
        }
    }

    //Restriction
    restriction(&r, &(*bTab)[niveauRel+1], mxc, myc, bornes, hc, h, mx-2, my-2);

    if(h>=hcMaxGrossier/2){
        //Résolution à la grille la plus grossière
                if (solve_umfpack((*nTab)[niveauRel+1], (*iaTab)[niveauRel+1], (*jaTab)[niveauRel+1], (*aTab)[niveauRel+1], (*bTab)[niveauRel+1], (*uTab)[niveauRel+1]))
            return 1;
    }
    else{
        //On déscend encore d'un niveau supplémentaire
        multigrille(uTab, nTab, 2*h, aTab, iaTab, jaTab, bTab, bornes, Lx, Ly, niveau, hcMaxGrossier, nbreItMG, rNorm);
    }

    //Prolongation
    prolongation(&(*uTab)[niveauRel], &(*uTab)[niveauRel+1], mxc, myc, bornes, h, hc, mx-2, my-2);
    
    //Calcul du résidu après la correction
    resVect(&(*uTab)[niveauRel], &(*aTab)[niveauRel], &(*iaTab)[niveauRel], &(*jaTab)[niveauRel], &(*bTab)[niveauRel], (*nTab)[niveauRel], &r);
    
    //Afficher résidu
    if(niveauRel==0){
        for(i=0; i<(*nTab)[0]; i++){
            rBis[i]=r[i]/sqrt(tot1);
        }
        if((nbreItMG==0 || nbreItMG==1) && afficherRes){
            printf("\n\nResidu avant post-smoothing \n\n");
            tableau = malloc(mx*my*sizeof(double));
            creaTableau(rBis, tableau, mx, my, bornes, h, 0.0);
            if (plot(rBis, mx, my, h, tableau, "dataResidu.txt", "fonctionResidu.txt"))
                return 1;
        }
    }
    
    //Post-smoothing
    for(i=0;i<4;i++){
        GS(&(*iaTab)[niveauRel], &(*jaTab)[niveauRel], &(*aTab)[niveauRel], &(*bTab)[niveauRel], &r, 0, &(*uTab)[niveauRel], (*nTab)[niveauRel]);
    }
    
    //Afficher résidu
    if(niveauRel==0){
        for(i=0; i<(*nTab)[0]; i++){
            rBis[i]=r[i]/sqrt(tot1);
        }
        if((nbreItMG==0 || nbreItMG==1) && afficherRes){
            printf("\n\nResidu après post-smoothing \n\n");
            tableau = malloc(mx*my*sizeof(double));
            creaTableau(rBis, tableau, mx, my, bornes, h, 0.0);
            if (plot(rBis, mx, my, h, tableau, "dataResidu.txt", "fonctionResidu.txt"))
                return 1;
        }
    }

    free(r);
    
    return 0;
}

