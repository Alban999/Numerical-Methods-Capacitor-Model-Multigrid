#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "prolongation.h"


int prolongation(double **u, double **uc, int mxc, int myc, double bornes[4], double h, double hc, int nx, int ny)
/*
 
 But
 ===
 Prolonger la grille grossière sur la grille fine
 
 Arguments
 =========
 u      (output) - pointeur vers le tableau 'u' (residu)
 uc     (input)  - pointeur vers le tableau 'rc' (residu restreint)
 mxc    (input)  - nombre de points dans la direction 1x sur la grille
 myc    (input)  - nombre de points dans la direction 1y sur la grille
 bornes (input)  - les bornes entre lesquelles est défini le condensateur
 h      (input)  - le pas de discrétisation de la grille fine
 hc     (input)  - le pas de discrétisation de la grille grossière
 nx     (input)  - le nombre d'éléments à l'intérieur du condensateur dans la directions 1x
 ny     (input)  - le nombre d'éléments à l'intérieur du condensateur dans la directions 1x
 
 */
{
    
    /* DECLARATION DES VARIABLES */
    
    int  nPNc = (int)((bornes[3]-bornes[2])/hc)+1; //Nombre de points de Neumann
    int xMinBefore = bornes[2]/h-1;
    int i, ixc = 0, iyc = 0, nxc = mxc-2, nyc = myc-2, indc = 0, ind = 0, ix = 0, iy = 0;
    
    /* PROLONGATION */
    
    for(iyc = 0; iyc < nyc; iyc++){
        for(ixc = 0; ixc < nxc; ixc++){
            indc = ixc + nxc * iyc;
            ix = 2*ixc+1;
            iy = 2*iyc+1;
            ind = ix + iy*nx;
            //SW
            (*u)[ind-nx-1] += 0.25*(*uc)[indc];
            //S
            (*u)[ind-nx] += 0.5*(*uc)[indc];
            //SE
            (*u)[ind-nx+1] += 0.25*(*uc)[indc];
            //W
            (*u)[ind-1] += 0.5*(*uc)[indc];
            //Sur le point
            (*u)[ind] += (*uc)[indc];
            //E
            (*u)[ind+1] += 0.5*(*uc)[indc];
            //NW
            (*u)[ind+nx-1] += 0.25*(*uc)[indc];
            //N
            (*u)[ind+nx] += 0.5*(*uc)[indc];
            //NE
            (*u)[ind+nx+1] += 0.25*(*uc)[indc];

        }
    }
    
    //Neumann
    for(ixc = 0; ixc < nPNc; ixc++){
        indc = ixc + nxc * nyc;
        ix = 2*ixc;
        ind = ix + ny*nx;

        //SW
        (*u)[ind-nx+xMinBefore-1] += 0.25*(*uc)[indc];
        //S
        (*u)[ind-nx+xMinBefore] += 0.5*(*uc)[indc];
        //SE
        (*u)[ind-nx+xMinBefore+1] += 0.25*(*uc)[indc];

        //Sur le point
        (*u)[ind] += (*uc)[indc];
        
        if(ixc == 0){
            //E
            (*u)[ind+1] += 0.5*(*uc)[indc];
        }
        else if(ixc == nPNc-1){
            //W
            (*u)[ind-1] += 0.5*(*uc)[indc];
        }
        else{
            //E
            (*u)[ind+1] += 0.5*(*uc)[indc];

            //W
            (*u)[ind-1] += 0.5*(*uc)[indc];
        }
    }
    
    
    
    return 0;
}

