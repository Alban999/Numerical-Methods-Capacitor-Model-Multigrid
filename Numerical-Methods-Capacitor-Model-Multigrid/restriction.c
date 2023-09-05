#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "restriction.h"

int restriction(double **r, double **rc, int mxc, int myc, double bornes[4], double hc, double h, int nx, int ny)
/*
 
 But
 ===
 Effectuer la restriction sur la grille fine
 
 Arguments
 =========
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
    
    int  nPNc = (int)((bornes[3]-bornes[2])/hc+1); //Nombre de points de Neumann
    int xMinBefore = bornes[2]/h-1;
    int ixc = 0, iyc = 0, nxc = mxc-2, nyc = myc-2, indc = 0, ind = 0, ix = 0, iy = 0;
    
    /* RESTRICTION */
    
    for(iyc = 0; iyc < nyc; iyc++){
        for(ixc = 0; ixc < nxc; ixc++){
            indc = ixc + nxc * iyc;
            ix = 2*ixc+1;
            iy = 2*iyc+1;
            ind = ix + iy*nx;
            
            (*rc)[indc] = 0.0625*(*r)[ind-nx-1] + 0.125*(*r)[ind-nx] + 0.0625*(*r)[ind-nx+1] + 0.125*(*r)[ind-1] + 0.25*(*r)[ind] + 0.125*(*r)[ind+1] + 0.0625*(*r)[ind+nx-1] + 0.125*(*r)[ind+nx] + 0.0625*(*r)[ind+nx+1];
        }
    }

    
    //Neumann
    ind = ny*nx;
    indc = nxc * nyc;
    (*rc)[indc] = 0.0625*(*r)[ind-nx+xMinBefore-1]+0.125*(*r)[ind-nx+xMinBefore]+0.0625*(*r)[ind-nx+xMinBefore+1]+0.25*(*r)[ind]+0.125*(*r)[ind+1];
    
    for(ixc = 1; ixc < nPNc-1; ixc++){
        indc = ixc + nxc * nyc;
        ix = 2*ixc;
        ind = ix + ny*nx;

        (*rc)[indc] = 0.0625*(*r)[ind-nx+xMinBefore-1]+0.125*(*r)[ind-nx+xMinBefore]+0.0625*(*r)[ind-nx+xMinBefore+1]+0.125*(*r)[ind-1]+0.25*(*r)[ind]+0.125*(*r)[ind+1];
    }
    ind = 2*(nPNc-1)+ny*nx;
    indc = nPNc-1 + nxc * nyc;
    (*rc)[indc] = 0.0625*(*r)[ind-nx+xMinBefore-1]+0.125*(*r)[ind-nx+xMinBefore]+0.0625*(*r)[ind-nx+xMinBefore+1]+0.125*(*r)[ind-1]+0.25*(*r)[ind];
    
    
    return 0;
}

