#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prob.h"

int prob(int mx, int my, double Lx, double Ly, int *n, int **ia, int **ja, double **a, double **b, double h, double bornes[4], double eps, double Q, int initializeRho, double (*rho)(), double uB)
/*
 
 But
 ===
 Générer le système linéaire n x n
 
 Au = b
 
 qui correspond à la discrétisation sur une grille cartésienne
 régulière m x m de l'équation de Poisson à deux dimensions
 
 d    d        d    d
 - == ( == u ) - == ( == u )  = 0     sur [0,1] x [0,1]
 dx   dx       dy   dy
 
 avec les conditions aux limites de Dirichlet
 
 u = 0  sur (x,0)                 , avec 0 <=  y  <= 1 ,
 u = 1  sur (0,y), (1,y) et (x,1) , avec 0 <= x,y <= 1 .
 
 La numérotation des inconnues est lexicographique, la direction x étant
 parcourue avant celle de y. La matrice A est retournée dans le format
 CRS qui est défini via trois tableaux : 'ia', 'ja' et 'a'.
 
 Arguments
 =========
 mx (input) - nombre de points dans la direction 1x sur la grille
 my (input) - nombre de points dans la direction 1y sur la grille
 (les valeurs de mx et my inférieures à 2 ne sont pas valides)
 Lx (input) - longeur du condensateur dans la direction 1x
 Ly (input) - longeur du condensateur dans la direction 1y
 n  (output) - pointeur vers le nombre d'inconnus dans le système
 ia (output) - pointeur vers le tableau 'ia' de la matrice A
 ja (output) - pointeur vers le tableau 'ja' de la matrice A
 a  (output) - pointeur vers le tableau 'a' de la matrice A
 b  (output) - pointeur vers le tableau 'b'
 h  (input) - le pas de discrétisation
 bornes (input) - les bornes entre lesquelles est défini le condensateur
 eps (input) - valeur de epsilon (la permittivité)
 Q  (input) - la charge du condensateur
 initializeRho - rendre rho nul (0 = Oui, 1 = Non)
 (*rho)(double, double, double, double, double, int) (input) - fonction qui renvoie la valeur de rho en fonction de où on se situe sur le condensateur
 uB (input) - Potentiel sur le bord du bas du condensateur compris entre les bornes données
 
 */
{
    
    /* DECLARATION DES VARIABLES */
    
    int  nnz, ix, iy, ind, nx, ny, nPN, xMinPN; // xMinPN : points avant les points nPN
    double invh2;
    
    nx = mx - 2;
    ny = my - 2;//CB ne sont pas pris en compte
    
    nPN = (int)((bornes[3]-bornes[2])/h+1); //Nombre de points de Neumann
    xMinPN = (int)(bornes[2]/h-1);

    
    //Verification des dimensions
    if(nx < 0) {
        printf("\n ERREUR : mx = %d n'est pas une valeur valide\n\n",mx);
        return 1;
    }
    else if(ny < 0) {
        printf("\n ERREUR : my = %d n'est pas une valeur valide\n\n",my);
        return 1;
    }
    
    invh2 = 1/(h*h); //h^-2 pour Ly=1
    *n  = nx * ny + nPN; //nombre d'inconnues
    nnz = 5 * nx * ny - 2 * nx - 2 * ny + 5 * nPN - 2; //nombre d'éléments non nuls

  
    //Allocation des tableaux
    *ia  = malloc((*n + 1) * sizeof(int));
    *ja  = malloc(nnz * sizeof(int));
    *a   = malloc(nnz * sizeof(double));
    *b   = malloc(*n * sizeof(double));
    
    //Vérification des allocations
    
    if (*ia == NULL || *ja == NULL || *a == NULL || *b == NULL ) {
        printf("\n ERREUR : pas assez de mémoire pour générer le système\n\n");
        return 1;
    }
    
    
    /* REMPLISSAGE DE LA MATRICE */
    
    ind = 0; /* au cas où nx = 0 */
    nnz = 0;
    
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {
            
            //Numéro de l'équation
            ind = ix + nx * iy;
            
            //Marquer le début de la ligne suivante dans le tableau 'ia'
            (*ia)[ind] = nnz;
            
            //Calculer le membre de droite
            (*b)[ind] = rho((ix+1)*h, (iy+1)*h, Lx, Ly, Q, initializeRho)/eps;
            //Remplissage de la ligne : voisin sud
            if (iy > 0)  {
                (*a)[nnz] = -invh2; //pour D=1
                (*ja)[nnz] = ind - nx;
                nnz++;
            }
            else if(ix>=bornes[0]/h-1 && ix<=bornes[1]/h-1)
                (*b)[ind] += uB*invh2;
            
            //Remplissage de la ligne : voisin ouest
            if (ix > 0)  {
                (*a)[nnz] = -invh2; //pour D=1
                (*ja)[nnz] = ind - 1;
                nnz++;
            } 
            
            //Remplissage de la ligne : élément diagonal
            (*a)[nnz] = 4.0*invh2; //pour D=1
            (*ja)[nnz] = ind;
            nnz++;
            
            //Remplissage de la ligne : voisin est
            if (ix < nx - 1) {
                (*a)[nnz] = -invh2; //pour D=1
                (*ja)[nnz] = ind + 1;
                nnz++;
            }
            
            //Remplissage de la ligne : voisin nord
            if (iy < ny - 1) {
                (*a)[nnz] = -invh2; //pour D=1
                (*ja)[nnz] = ind + nx;
                nnz++;
            }
            else if (ix>= bornes[2]/h-1 && ix<= bornes[3]/h-1){
                (*a)[nnz] = -invh2; //pour D=1
                (*ja)[nnz] = ind - xMinPN + nx;
                nnz++;
            }
            
        }
    }

    
    for (ix = 0; ix < nPN; ix++) {
        
        //Numéro de l'équation
        ind = ix + nx * ny;
        
        //Marquer le début de la ligne suivante dans le tableau 'ia'
        (*ia)[ind] = nnz;
        //Calculer le membre de droite
        (*b)[ind] = rho((ix+1)*h+bornes[2], Ly, Lx, Ly, Q, initializeRho)/(2*eps);
        
        
        //Remplissage de la ligne : voisin sud
        (*a)[nnz] = -invh2; //pour D=1
        (*ja)[nnz] = xMinPN + ind - nx;
        
        nnz++;
        
        //Remplissage de la ligne : voisin ouest
        if (ix > 0)  {
            (*a)[nnz] = -0.5*invh2; //pour D=1
            (*ja)[nnz] = ind - 1;
            nnz++;
        }
        //Remplissage de la ligne : élément diagonal
        (*a)[nnz] = 2.0*invh2; //pour D=1
        (*ja)[nnz] = ind;
        nnz++;
        //Remplissage de la ligne : voisin est
        if (ix < nPN - 1) {
            (*a)[nnz] = -0.5*invh2; //pour D=1
            (*ja)[nnz] = ind + 1;
            nnz++;
        }
        
    }

    
    //dernier élément du tableau 'ia'
    (*ia)[ind + 1] = nnz;
    
    return 0;
}

