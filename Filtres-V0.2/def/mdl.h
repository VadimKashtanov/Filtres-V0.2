#pragma once

#include "marchee.h"

//	Instructions:
//	0|	n-filtre-ema-prixs
//	1|	n-filtre
//	2|	n-neurone

#define POIDS_NEU(n) (2*n+1)
#define CONSTS_FLTR(n) (n)

#define LOCDS_NEU(n) (1+n)
#define LOCDS_FLTR(n) (6+n+(n-1))

typedef struct {
	//	Structure du systeme
	uint couches, N, C;
	uint * type;	// dans {0,1,2}, la 0-eme etant toujours de 0
	uint * y;	//	combien de sorties
	uint * n;	//	combien d'entrees prend chaque y

	uint max_n;

	//	Couches neuronalles
	uint *** neu_vers;	//[couche][neurone][connection] pris dans {0..y[-1]}
	uint ** fltr_depart;	//[couche][filtre] pris dans {0..(y[-1] - n)}

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
	uint * couche_type, uint * n, uint * y,
	uint *** neu_vers,	//[i] = 0 si i-eme est non neuronal
	uint ** fltr_depart,	//[i] = 0 si ieme est non filtrique
	// 0 -eme couche
	uint * _intervalles,
	uint * _ema);
void liberer_mdl(Mdl_t * mdl);

void verifier_derivee(Mdl_t * mdl);

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
float score(Mdl_t * mdl);
//=somme(gains); w -= f'(x) * alpha; alpha = 1/moy(ddf(x))