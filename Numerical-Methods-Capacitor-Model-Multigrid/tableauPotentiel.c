#include <stdio.h>
#include <stdlib.h>
#include "tableauPotentiel.h"

int creaTableau(double *x, double *tableauPot, int mx, int my, double bornes[4], double h, double uB)
/*
 
 But
 ===
 Créer un tableau qui contient toutes les valeurs du potentiel à chaque de la grille.
 
 Arguments
 =========
 x  (input) - le vecteur solution
 tableauPotentiel (input) - tableau composé des différentes valeurs du potentiel en fonction de la position dans le condensateur
 mx (input) - nombre de points dans la direction 1x sur la grille
 my (input) - nombre de points dans la direction 1y sur la grille
 bornes (input) - les bornes entre lesquelles est défini le condensateur
 h  (input) - le pas de discrétisation
 uB (input) - Potentiel sur le bord du bas du condensateur compris entre les bornes données
 
 */
{
    int i, j, tot=0, tt=0;
    
    /* CREATION DU TABLEAU DE POTENTIEL ENTIER SUR TOUT LE CONDENSATEUR */
    
    for(i=0;i<mx;i++){
        if(i*h>=bornes[0] && i*h<=bornes[1])
            tableauPot[tot] = uB;
        else
            tableauPot[tot] = 0.0;
        tot++;
    }

    for (i=0; i<my-2; i++) {
        for(j=0;j<mx;j++)
        {
            if(j==0||j==mx-1)
                tableauPot[tot] = 0.0;
            else{
                tableauPot[tot] = (x)[tt];
                tt++;
            }
            tot++;
        }
    }
    
    
    for(i=0;i<mx;i++){
        if(i*h>=bornes[2] && i*h<=bornes[3]){
            tableauPot[tot] = (x)[tt];
            tt++;
        }
        else
            tableauPot[tot] = 0.0;
        tot++;
    }
    
    return 0;

}
    
