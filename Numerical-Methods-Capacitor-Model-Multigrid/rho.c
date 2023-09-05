#include <stdio.h>
#include <math.h>
#include "rho.h"

double rho(double x, double y, double Lx, double Ly, double Q, int initialize)
/*
 
 But
 ===
 Calculer la densité électrique au sein du condensateur,
 sachant que la distribution de charge correspond à un
 disque uniformément chargé, centré au milieu de la
 surface Ω du condensateur et de rayon Ly/4.
 
 Formule d'un disque :
 (x-xc)^2 * (y-yc)^2 <= R^2
 
 Arguments
 =========
 x  (input) - position sur l'axe des x (en mètres)
 y  (input) - position sur l'axe des y (en mètres)
 Lx (input) - longeur du condensateur dans la direction 1x
 Ly (input) - longeur du condensateur dans la direction 1y
 Q  (input) - la charge du condensateur
 initialize - rendre rho nul (0 = Oui, 1 = Non)

 */
{
    /* DECLARATION DES VARIABLES */
    
    double resultat, r, val, xc, yc, rhoResult;
    
    if(initialize){
        //rayon du disque
        r = Ly/4;
        
        //Centre du disque
        xc = Lx/2;
        yc = Ly/2;
        
        //Définir valeur de rho par unité de valeur selon z
        val = Q/(M_PI*r*r);
        
        
        /* CALCUL DE RHO */
        
        resultat = (x-xc)*(x-xc) + (y-yc)*(y-yc);

        if(resultat <= r*r)
            rhoResult = val;
        else
            rhoResult = 0.0;
    }
    else{
        rhoResult = 0.0;
    }
    return rhoResult;
}
