#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "residuEvolution.h"

int residuEvolution(int nbreIter, double **rNorm, char *nameFile, char *titre)
/*
 
 But
 ===
 Représenter l'évolution du résidu en fonction du nombre d'itérations.
 
 Arguments
 =========
 nbreIter (input) - nombre d'itérations effectué
 rNorm    (input) - tableau contenant les différents résidus pour chaque itération
 nameFile (input) - nom du fichier contenant les données
 titre    (input) - titre du graphe
 
 */
{
    
    //Ouverture du fichier de données
    FILE *fichierDonnees;
    fichierDonnees = fopen(nameFile, "w");
    for(int i = 0; i<nbreIter; i++){
       
        fprintf(fichierDonnees, "%d %f\n", i+1, log((*rNorm)[i]));
        
    }
    
    fclose(fichierDonnees);

    //Création du fichier de commandes
    FILE *cmd;
    cmd = popen("gnuplot -persistent", "w");
    
    fprintf(cmd, "load 'fonctionResiduEnFctDuNbreIter.txt'\n");
    fprintf(cmd, "set title \"Etude du résidu en fonction du nombre d'itérations : %s\"\n", titre);
    fprintf(cmd, "set xrange[%d:%d]\n", 0, nbreIter);
    fprintf(cmd, "set yrange[%f:%f]\n", log((*rNorm)[nbreIter-1]), log((*rNorm)[0]));
    fprintf(cmd, "plot 'dataResiduEnFctDuNbreIter.txt' using 1:2 with linespoints title 'TWO GRID'\npause -1");
    
    //Vérification de la création des fichiers
    if (cmd == NULL || fichierDonnees == NULL)
    {
        printf("Erreur d'ouverture de fichier pour la création de la fonction du résidu en fonction du pas de discrétisation.\n");
        return 1;
    }

    pclose(cmd);
    

    return 0;
}
