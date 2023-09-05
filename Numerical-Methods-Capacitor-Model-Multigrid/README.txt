Alban Dietrich, Master-Physique, n°ID : 000440726, 12/12/2019

Dans le main :

- Sélection du nombre de niveau avec la variable ‘niv’

- Choix d’afficher les graphes des résidus (visualisation) à l’aide de la variable ‘afficherRes’

- Définition du nombre d'itérations des cycles : 1. ‘nbreIterTG’ -> Two-Grid
						 2. ‘nbreIterMG’ -> Multigrid
						 3. ‘nbreIterGC’ -> Gradient conjugué 


- Choix de l’exercice à l’aide des variables : 1. ‘Exercice 1’ -> Neumann
					       2. ‘Exercice 2’ -> Two-Grid
					       3. ‘Exercice 3’ -> Multigrid
					       4. ‘Exercice 4’ -> Gradient conjugué


Remarque: Je n’ai pas réussi à faire fonctionner correctement mon code de gradient conjugué avec le Multigrid comme préconditionneur. J’ai laissé le code dans le projet si vous voulez voir ce qu’il affiche. 
Par contre j’ai réussi à faire un code pour calculer le gradient conjugué sans préconditionneur. Ce fichier est ‘CG_Sans_Preconditionneur.c’, je le mets dans le dossier du projet sans le relier au Makefile si vous voulez y jeter un coup d’oeil.