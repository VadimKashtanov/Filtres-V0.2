#include "mdl.h"

static uint* cpyuint(uint * arr, uint len) {
	uint * ret = malloc(sizeof(uint) * len);
	memcpy(ret, arr, sizeof(uint) * len);
	return ret;
}

static float* allouer_flotants(uint nb) {
	return malloc(sizeof(float) * nb);
}

Mdl_t * cree_model(
	uint couches,
	uint * couche_type, uint * couche_n, uint * couche_y,
	uint *** couche_neurone_conn,	//[i] = 0 si i-eme est non neuronal
	uint ** couche_filtre_depart,	//[i] = 0 si ieme est non filtrique
	// 0 -eme couche
	uint * _intervalles,
	uint * _ema)
{
	assert(couches > 0);
	assert(couche_type[0] == 0);
	assert(couche_y[couches-1] == 1);
	//
	Mdl_t * mdl = malloc(sizeof(Mdl_t));
	//
	uint N = couches;
	mdl->N = N;
	mdl->couches = couches;
	mdl->couche_type = cpyuint(couche_type, N);
	mdl->couche_n = cpyuint(couche_n, N);
	mdl->couche_y = cpyuint(couche_y, N);
	//
	mdl->intervalles = cpyuint(_intervalles, couche_n[0]);
	mdl->ema = cpyuint(_ema, couche_n[0]);
	//
	mdl->locds = 0;
	mdl->vars = 0;
	mdl->locd_depart = malloc(sizeof(uint) * N);
	mdl->y_depart = malloc(sizeof(uint) * N);
	mdl->poids = 0;
	mdl->constantes = 0;
	mdl->conste_depart = malloc(sizeof(uint) * N);
	mdl->poid_depart = malloc((sizeof)uint * N);
	//
	mdl->couche_neurone_conn = calloc(N, sizeof(uint**));
	mdl->couche_filtre_depart = calloc(N, sizeof(uint*));
	for (uint i=1; i < N; i++) {
		if (couche_type[i] == 2) {
			//	=== Couches Neurones ===
			mdl->couche_neurone_conn[i] = malloc(sizeof(uint*) * couche_y[i]);
			for (uint j=0; j < couche_y[i]; j++)
				mdl->couche_neurone_conn[i][j] = cpyuint(couche_neurone_conn[i][j], couche_n[i])
			//
			mdl->poid_depart[i] = mdl->poids;
			mdl->conste_depart[i] = mdl->constantes;
			mdl->locd_depart[i] = mdl->locds;
			//
			mdl->poids += couche_y[i] * (2*couche_n[i] + 1);
			mdl->constantes += 0;
			mdl->locds += 6+couche_n[i]+couche_n[i]-1;
		} else {
			// === Couches Filtriques ===
			mdl->couche_filtre_depart[i] = cpyuint(couche_filtre_depart[i], couche_y[i]);
			//
			mdl->poid_depart[i] = mdl->poids;
			mdl->conste_depart[i] = mdl->constantes;
			//
			mdl->poids += 0;
			mdl->constantes += couche_y[i]*couche_n[i];
			mdl->locds += 1;
		}
		//
		mdl->y_depart[i] = mdl->vars;
		mdl->vars += couche_y[i];
	}
	//
	mdl->var = allouer_flotants(mdl->vars);
	mdl->d_var = allouer_flotants(mdl->vars);
	mdl->d_poid = allouer_flotants(mdl->poids);
	mdl->constante = allouer_flotants(mdl->constantes);
	mdl->poid = allouer_flotants(mdl->poids);
	mdl->locd = allouer_flotants(mdl->locds);
	//
	return mdl;
}

void liberer_mdl(Mdl_t * mdl) {
	uint N = mdl->N;
	//
	for (uint i=1; i < N; i++) {
		if (couche_type[i] == 2) {
			for (uint j=0; j < couche_y[i]; j++)
				free(mdl->couche_neurone_conn[i][j]);
			free(mdl->couche_neurone_conn[i]);
		} else {
			free(mdl->couche_filtre_depart[i]);

		}
	}
	free(mdl->couche_neurone_conn);
	free(mdl->couche_filtre_depart);
	//
	free(mdl->couche_n);
	free(mdl->couche_y);
	free(mdl->couche_type);
	//
	free(mdl->intervalles);
	free(mdl->ema);
	free(mdl->constante);
	free(mdl->poid);
	//
	free(mdl->var);
	free(mdl->d_var);
	free(mdl->d_poid);
	free(mdl->conste_depart);
	free(mdl->poid_depart);
	free(mdl->y_depart);
	free(mdl->locd);
	free(mdl->locd_depart);
	//
	free(mdl);
};