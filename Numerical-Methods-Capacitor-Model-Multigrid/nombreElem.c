#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "nombreElem.h"

int nombreElem(double h, double Lx, double Ly, double bornes[4])
/*
 
 But
 ===
 Renvoyer le nombre d'éléments de la grille en fonction du pas de discrétisation
 
 Arguments
 =========
 h       (input) - le pas de discrétisation
 Lx      (input) - longeur du condensateur dans la direction 1x
 Ly      (input) - longeur du condensateur dans la direction 1y
 bornes  (input) - les bornes entre lesquelles est défini le condensateur
 
 */
{
    int my = Ly/h+1;
    int mx = Lx/h+1;
    int nx = mx - 2;
    int ny = my - 2;
    int nPN = (int)((bornes[3]-bornes[2])/h+1);
    
    return nx * ny + nPN;
    
}
