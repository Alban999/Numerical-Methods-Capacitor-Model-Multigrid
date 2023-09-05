#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "prob.h"
#include "umfpk.h"
#include "time.h"
#include "residu.h"
#include "tableauPotentiel.h"
#include "rho.h"
#include "GS.h"
#include "restriction.h"
#include "prolongation.h"
#include "multigrille.h"
#include "nombreElem.h"
#include "resVect.h"
#include "residuEvolution.h"
#include "residu.h"
#include "conjGrad.h"

int main(int argc, char *argv[])
/*
 
 But
 ===
 Gérer l'ensemble des fonctionnalités du projet
 
 */
{
    
    /* DECLARER LES VARIABLES */
    
    //Définition du nombre de niveaux
    int niv = 5;
    
    //Dimensions du condensateur
    double Lx = 0.01, Ly = 0.004;
    double hcMaxGrossier = 0.001;
    double h = hcMaxGrossier/pow(2, niv);
    
    //Définition du pas grossier pour la méthode TWO-GRID
    double hc  = h*2;

    //Faire en sorte que mx et my soient impairs
    int my = Ly/h+1, myc = Ly/hc+1, mx = Lx/h+1, mxc = Lx/hc+1;
    
    //Type de condensateur
    double epsVal = 8.5;
    double eps0 = 8.854*1e-12;
    double eps = epsVal*eps0;
    double Q = -0.2e-7;
    int rhoInitialize = 0; //Si = 0 -> rho = 0 sinon rho != 0

    //Les bornes du condensateur
    //Faire en sorte que les bornes supérieurs ne soient pas sur le bord du condensateur pour les CB de Neumann
    //Pour les CB de Neuman,  on oblige à ce que les extrémités de uc tombe sur les extremités de u
    double bornes[4] = {0.002, 0.008, 0.002, 0.004};
    
    //Potentiel sur le bord bas du condensateur
    double uB = -0.5;
    
    //Format CSR
    int *ia, *ja, *iac, *jac;
    double *a, *b, *x, *r, *rBis, *u, *uc, *ac, *bc, *rc, *d, *u0;
    
    //Les temps calculés avant et après la résolution du solveur
    double t1, t2;
    
    //Tableau à afficher graphiquement
    double *tableau;
    
    //Autres variables
    int n, nc, i, niveau = -1, p = 0;
    double tot1 = 0;
    
    //Définition du nombre d'itérations des cycles
    int nbreIterTG = 50, nbreIterMG = 50, nbreIterGC = 100;
    
    //Afficher les graphes de résidus : oui ou non
    int afficherRes = 0;
    
    //Temps de relaxation
    double tau = 1.0;
    
    //Vecteur norme residu
    double *rNormTG, *rNormMG;
    
    //Création des tableaux de pointeurs pour la fonction multigrille
    double **uTab, **aTab, **bTab;
    int *nTab, **iaTab, **jaTab;
    
    //Nombre de niveaus + nombre d'éléments par niveau
    int niveauTot = log(hcMaxGrossier/h)/log(2)+1, nElem = 0;
    
    //Déclaration des tableaux de pointeurs pour le multigrille
    uTab = malloc(niveauTot * sizeof(*uTab));
    aTab = malloc(niveauTot * sizeof(*aTab));
    bTab = malloc(niveauTot * sizeof(*bTab));
    iaTab = malloc(niveauTot * sizeof(*iaTab));
    jaTab = malloc(niveauTot * sizeof(*jaTab));
    nTab = malloc(niveauTot * sizeof(int));
    
    if(uTab==NULL || aTab==NULL || bTab==NULL || iaTab==NULL || jaTab==NULL || nTab==NULL){
        printf("\nErreur : problème d'allocation\n");
        return 1;
    }
    
    for(i = 0; i<niveauTot; i++){
        nElem = nombreElem(pow(2, i)*h, Lx, Ly, bornes);
        uTab[i] = malloc(nElem*sizeof(double));
        if(uTab[i]==NULL){
            printf("\nErreur : problème d'allocation\n");
            return 1;
        }
    }
    
    //Choix des exercices : oui ou non
    int Exercice1 = 0, Exercice2 = 1, Exercice3 = 0, Exercice5 = 0;
    
    
    /* GENERER LE PROBLEME */

    printf("\n|GENERATION DU PROBLEME|\n");

    if (prob(mx, my, Lx, Ly, &n, &ia, &ja, &a, &b, h, bornes, eps, Q, rhoInitialize, &rho, uB))
        return 1;

    printf("\nPROBLEM: ");
    printf("mx = %d my = %d  n = %d  nnz = %d\n", mx, my, n, ia[n]);


    if(Exercice1){
        
        /* SOLUTION AVEC CB DE NEUMANN */
        
        x = malloc(n * sizeof(double));
        
        if(x==NULL){
            printf("\nErreur : problème d'allocation\n");
            return 1;
        }
        
        //Calcul de la solution
        t1 = mytimer();
        if (solve_umfpack(n, ia, ja, a, b, x))
            return 1;
        t2 = mytimer();
        
        //Temps de résolution
        printf("\nTemps de solution UMFPACK (CPU): %5.1f sec\n",t2-t1);
        
        //Création du tableau de potentiel + génération du graphe
        tableau = malloc(mx*my*sizeof(double));
        creaTableau(x, tableau, mx, my, bornes, h, uB);
        if (plot(x, mx, my, h, tableau, "dataPotentiel.txt", "fonctionPotentiel.txt"))
            return 1;
        
        /* LIBERATION DE LA MEMOIRE*/
        
        free(x);
        
    }

    if(Exercice2){
        
        /* TWO GRID METHOD */
        
        u = malloc(n * sizeof(double));
        r = malloc(n * sizeof(double));
        rBis = malloc(n * sizeof(double));
        rNormTG = malloc(nbreIterTG * sizeof(double));
        
        for(i=0; i<n; i++){
            u[i] = 0.0;
            r[i] = b[i];
            tot1 += b[i]*b[i];
        }
        
        //Création de Ac
        if (prob(mxc, myc, Lx, Ly, &nc, &iac, &jac, &ac, &bc, hc, bornes, eps, Q, rhoInitialize, &rho, uB))
            return 1;
        
        rc = malloc(nc * sizeof(double));
        uc = malloc(nc * sizeof(double));

        for(p=0; p<nbreIterTG; p++){
            
            printf("\nItération %d\n", p);
            
            //Afficher résidu
            if((p==0 || p==1) && afficherRes){
                for(i=0; i<n; i++){
                    rBis[i]=r[i]/sqrt(tot1);
                }
                printf("\n\nResidu avant pré-smoothing \n\n");
                tableau = malloc(mx*my*sizeof(double));
                creaTableau(rBis, tableau, mx, my, bornes, h, 0.0);
                if (plot(rBis, mx, my, h, tableau, "dataResidu.txt", "fonctionResidu.txt"))
                    return 1;
            }
            
            //Pré-smoothing
            for(i=0; i<1; i++){
                GS(&ia, &ja, &a, &b, &r, 1, &u, n);
            }
            
            
            //Afficher résidu
            if((p==0 || p==1) && afficherRes){
                for(i=0; i<n; i++){
                    rBis[i]=r[i]/sqrt(tot1);
                }
                printf("\n\nResidu après pré-smoothing \n\n");
                tableau = malloc(mx*my*sizeof(double));
                creaTableau(rBis, tableau, mx, my, bornes, h, 0.0);
                if (plot(rBis, mx, my, h, tableau, "dataResidu.txt", "fonctionResidu.txt"))
                    return 1;
            }
            
            
            //Restriction
            restriction(&r, &rc, mxc, myc, bornes, hc, h, mx-2, my-2);

            //Calcul de la solution uc
            if (solve_umfpack(nc, iac, jac, ac, rc, uc))
                return 1;

            //Prolongation
            prolongation(&u, &uc, mxc, myc, bornes, h, hc, mx-2, my-2);

            //Calcul du résidu après la correction
            resVect(&u, &a, &ia, &ja, &b, n, &r);
            
            //Afficher résidu
            if((p==0 || p==1) && afficherRes){
                for(i=0; i<n; i++){
                    rBis[i]=r[i]/sqrt(tot1);
                }
                printf("\n\nResidu avant post-smoothing \n\n");
                tableau = malloc(mx*my*sizeof(double));
                creaTableau(rBis, tableau, mx, my, bornes, h, 0.0);
                if (plot(rBis, mx, my, h, tableau, "dataResidu.txt", "fonctionResidu.txt"))
                    return 1;
            }
            
            //Post-smoothing
            for(i=0; i<1; i++){
                GS(&ia, &ja, &a, &b, &r, 0, &u, n);
            }
            
            //Afficher résidu
            if((p==0 || p==1) && afficherRes){
                for(i=0; i<n; i++){
                    rBis[i]=r[i]/sqrt(tot1);
                }
                printf("\n\nResidu après post-smoothing \n\n");
                tableau = malloc(mx*my*sizeof(double));
                creaTableau(rBis, tableau, mx, my, bornes, h, 0.0);
                if (plot(rBis, mx, my, h, tableau, "dataResidu.txt", "fonctionResidu.txt"))
                    return 1;
            }
            
            //Stockage du résidu à chaque itération
            rNormTG[p] = residu(&u, &a, &ia, &ja, &b, n)*1e16;
            
        }
        
        //Création du tableau de potentiel + génération du graphe de la solution
        tableau = malloc(mx*my*sizeof(double));
        creaTableau(u, tableau, mx, my, bornes, h, uB);
        if (plot(u, mx, my, h, tableau, "dataPotentiel.txt", "fonctionPotentiel.txt"))
            return 1;

        //Affichage de l'évolution du résidu en fonction du nombre d'itérations
        if(residuEvolution(nbreIterTG, &rNormTG, "dataResiduEnFctDuNbreIter.txt", "TWO GRID"))
            return 1;
        
        /* LIBERATION DE LA MEMOIRE */
        
        free(iac); free(jac); free(ac); free(bc); free(uc); free(rc); free(rBis);free(r);
        
    }


    if(Exercice3){
        
        /* MULTIGRILLES */
        
        //Création de tous les tableaux pour chaque niveau
        double hBis;
        int myBis, mxBis;
        
        for(int m = 0; m<niveauTot; m++){
            hBis = pow(2, m)*h;
            myBis = Ly/hBis+1;
            mxBis = Lx/hBis+1;
            
            if (prob(mxBis, myBis, Lx, Ly, &nTab[m], &iaTab[m], &jaTab[m], &aTab[m], &bTab[m], hBis, bornes, eps, Q, rhoInitialize, &rho, uB))
                return 1;
        }
        
        //Allocation du tableau des normes des résidus en fonction du nombre d'itérations
        rNormMG = malloc(nbreIterMG * sizeof(double));
        
        if(rNormMG==NULL){
            printf("\nErreur : problème d'allocation\n");
            return 1;
        }
        
        //Initialisation
        for(int i =0; i< niveauTot; i++){
            for (int j = 0; j< nTab[i]; j++){
                uTab[i][j] = 0.0;
            }
        }
        
        u0 = malloc(nTab[0]*sizeof(double));
        
        if(u0==NULL){
            printf("\nErreur : problème d'allocation\n");
            return 1;
        }
        
        for(int w = 0; w<nbreIterMG; w++){
            
            printf("\nItération %d\n", w);
            
            //Réinitialise le niveau
            niveau = -1;
            
            for(i=0; i<nTab[0]; i++){
                u0[i] = uTab[0][i];
            }
            
            multigrille(&uTab, &nTab, h, &aTab, &iaTab, &jaTab, &bTab, bornes, Lx, Ly, &niveau, hcMaxGrossier, w, &rNormMG);
        
            for(i=0; i<nTab[0]; i++){
                uTab[0][i] = u0[i] + tau*(uTab[0][i] -u0[i]);
            }
            
            for(int i = 1; i< niveauTot; i++){
                for (int j = 0; j< nTab[i]; j++){
                    uTab[i][j] = 0.0;
                    bTab[i][j] = 0.0;
                }
            }
            
            //Stockage de la norme du résidu à chaque itération
            rNormMG[w] = residu(&uTab[0], &aTab[0], &iaTab[0], &jaTab[0], &bTab[0], nTab[0])*1e16;
            
            printf("\nresidu %e\n", rNormMG[w]*1e-16);

        }
        
  

        
        //Création du tableau de potentiel + génération du graphe
        tableau = malloc(mx*my*sizeof(double));
        creaTableau(uTab[0], tableau, mx, my, bornes, h, uB);
        if (plot(uTab[0], mx, my, h, tableau, "dataPotentiel.txt", "fonctionPotentiel.txt"))
            return 1;

        //Affichage de l'évolution du résidu en fonction du nombre d'itérations
        if(residuEvolution(nbreIterMG, &rNormMG, "dataResiduEnFctDuNbreIter.txt", "MULTIGRID"))
            return 1;
        
        
    }
    
    if(Exercice5){
        
        /* GRADIENT CONJUGUE */
        
        //Création de tous les tableaux pour chaque niveau
        double hBis;
        int myBis, mxBis;
        
        for(int m = 0; m<niveauTot; m++){
            hBis = pow(2, m)*h;
            myBis = Ly/hBis+1;
            mxBis = Lx/hBis+1;
            
            if (prob(mxBis, myBis, Lx, Ly, &nTab[m], &iaTab[m], &jaTab[m], &aTab[m], &bTab[m], hBis, bornes, eps, Q, rhoInitialize, &rho, uB))
                return 1;
            
        }
        
        u = malloc(nTab[0]*sizeof(double));
        
        if(u==NULL){
            printf("\nErreur : problème d'allocation\n");
            return 1;
        }
        
        //Allocation du tableau des normes des résidus en fonction du nombre d'itérations
        rNormMG = malloc(nbreIterMG * sizeof(double));
        
        //Gradient conjugué
        conjGrad(&iaTab, &jaTab, &aTab, &bTab, &u, &uTab, &nTab, h, bornes, Lx, Ly, &niveau, hcMaxGrossier, i, &rNormMG, nbreIterGC);
        
        //Création du tableau de potentiel + génération du graphe
        tableau = malloc(mx*my*sizeof(double));
        creaTableau(u, tableau, mx, my, bornes, h, uB);
        if (plot(u, mx, my, h, tableau, "dataPotentiel.txt", "fonctionPotentiel.txt"))
            return 1;
        
        
    }
    
    /* LIBERATION DE LA MEMOIRE */
    
    free(ia); free(ja); free(a); free(b); free(u); free(tableau); free(nTab);

    if(Exercice3 || Exercice5){
        for(i = 0; i<niveauTot; i++){
            free(uTab[i]); free(iaTab[i]); free(jaTab[i]); free(aTab[i]); free(bTab[i]);
        }
    }
    
  
    
    return 0;
}



