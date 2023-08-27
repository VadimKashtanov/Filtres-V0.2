#pragma once

#include "marchee.h"

//	Instructions:
//	0|	n-filtre-ema-prixs
//	1|	n-filtre
//	2|	n-neurone

typedef struct {
	//	Structure du systeme
	uint couches, N;
	uint * couche_type;	// dans {0,1,2}, la 0-eme etant toujours de 0
	uint * couche_n;	//	combien d'entrees prend chaque y
	uint * couche_y;	//	combien de sorties

	//	Couches neuronalles
	uint *** couche_neurone_conn;	//[couche][neurone][connection] pris dans {0..y[-1]}
	uint ** couche_filtre_depart;	//[couche][filtre] pris dans {0..(y[-1] - n)}

	//	Pour la 0-eme couche
	uint * intervalles;	// dans intervalles[INTERVALLES] (voir marchee.h)
	uint * ema;			// dans {0, 3, 5, 10, 25, 50, 100, 250, 500, 1000} (voir marchee.h)

	//	Les poids mutables par selection ou pas optimisation analytique
	uint constantes, poids;
	float * constante;	//selection naturelle
	float * poid;		// -f'(x)
 
	//	Les espaces de calcule F(x) et retro-propagation F'(x)
	uint vars, locds;
	float * var;
	float * locd;
	float * d_var, * d_poid;

	//	A partire d'ou chaque couche prend ces constantes et ces poids
	uint * conste_depart, * poid_depart, * y_depart, * locd_depart;
} Mdl_t;

//	Allocation Memoire
Mdl_t * cree_mdl(
	uint couches,
	uint * couche_type, uint * couche_n, uint * couche_y,
	uint *** couche_neurone_conn,	//[i] = 0 si i-eme est non neuronal
	uint ** couche_filtre_depart,	//[i] = 0 si ieme est non filtrique
	// 0 -eme couche
	uint * _intervalles,
	uint * _ema);
void liberer_mdl(Mdl_t * mdl);

//	Disque Dur
void ecrire_mdl(Mdl_t * mdl, char * fichier);
Mdl_t * lire_mdl(char * fichier);

//	Etc
void plume_mdl(Mdl_t * mdl);

//	F(x) & F'(x)
float f(Mdl_t * mdl, uint depart);
void df(Mdl_t * mdl, uint depart, float erreur);

//	Fonction Objectif
float objectif_gain(Mdl_t * mdl, uint depart);
void d_objectif_gain(Mdl_t * mdl, uint depart, float obj_gain);

//	Optimisation & Gain
float score(Mdl_t * mdl);	//=somme(gains); w -= f'(x) * alpha; alpha = 1/moy(ddf(x))