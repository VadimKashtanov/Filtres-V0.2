#include "marchee.h"

float prixs[PRIXS] = {};
float ema[NB_DIFF_EMA][PRIXS] = {};

uint intervalles[INTERVALLES] = {
	1,1,1,1,1,
	2,3,4,5,6,
	10,15,25,50,75,
	100,200,300,500,1000
};

static uint __ema[NB_DIFF_EMA] = {
	0, 3, 5, 12, 26, 50, 100, 250, 500, 1000
};

float probas_choix_ema[NB_DIFF_EMA] = {0};

static float probas_ema[NB_DIFF_EMA] = {	//Assez arbitraire, le but est de faire donner plus d'analyse de em0 que les autres
	2,.5,.3, .3, .3, .3,  .3,  .3,  .5,   .5
};

void charger_les_prixs() {
	//	Au passage on calcule aussi ```probas_choix_ema```
	float _sum = 0;
	for (uint i=0; i < NB_DIFF_EMA; i++) _sum += probas_ema[i];
	//
	probas_choix_ema[0] = 0.0;
	for (uint i=1; i < NB_DIFF_EMA; i++) probas_choix_ema[i] = probas_choix_ema[i-1] + probas_ema[i-1]/_sum;
	

	//	Chargement des Prixs et calcule des EMA, Moyennes Mouvantes
	FILE * fp = fopen("prixs/prixs.bin", "rb");
	assert(fp != 0);
	uint __PRIXS;
	fread(&__PRIXS, sizeof(uint), 1, fp);
	assert(__PRIXS == PRIXS);
	fread(prixs, sizeof(float), PRIXS, fp);
	fclose(fp);
	//
	float k[NB_DIFF_EMA];
	float _k[NB_DIFF_EMA];
	for (uint i=0; i < NB_DIFF_EMA; i++) {
		k[i] = 1.0/(1.0 + (float)__ema[i]);
		_k[i] = 1.0 - k[i];
		ema[i][0] = prixs[0];
	}
	//
	for (uint i=1; i < PRIXS; i++) {
		for (uint j=0; j < NB_DIFF_EMA; j++) {
			ema[j][i] = prixs[i]*k[j] + ema[j][i-1]*_k[j];
		};
	};
};