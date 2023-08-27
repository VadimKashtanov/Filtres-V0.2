#include "mdl.h"

static uint* cpyuint(uint * arr, uint len) {
	uint * ret = malloc(sizeof(uint) * len);
	memcpy(ret, arr, sizeof(uint) * len);
	return ret;
}

static float* allouer_flotants(uint nb) {
	return malloc(sizeof(float) * nb);
}

Mdl_t * cree_mdl(
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
	assert(couche_n[0] <= MAX_N_PREMIERE_COUCHE);
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
	mdl->intervalles = cpyuint(_intervalles, couche_y[0]);
	mdl->ema = cpyuint(_ema, couche_y[0]);
	//
	for (uint i=0; i < couche_y[0]; i++) {
		assert(mdl->intervalles[i] <= intervalles[INTERVALLES-1]);
		assert(mdl->ema[i] <= NB_DIFF_EMA);
	};
	//
	mdl->locds = 0;
	mdl->vars = 0;
	mdl->locd_depart = malloc(sizeof(uint) * N);
	mdl->y_depart = malloc(sizeof(uint) * N);
	mdl->poids = 0;
	mdl->constantes = 0;
	mdl->conste_depart = malloc(sizeof(uint) * N);
	mdl->poid_depart = malloc(sizeof(uint) * N);
	//
	mdl->couche_neurone_conn = calloc(N, sizeof(uint**));
	mdl->couche_filtre_depart = calloc(N, sizeof(uint*));
	//
	//	==== Couche 0 =====
	mdl->poid_depart[0] = 0;
	mdl->conste_depart[0] = 0;
	mdl->locd_depart[0] = 0;
	mdl->y_depart[0] = 0;
	//
	mdl->vars = couche_y[0];
	mdl->poids = 0;
	mdl->constantes = couche_n[0]*couche_y[0];
	mdl->locds = couche_y[0]*(6+couche_n[0]+couche_n[0]-1);
	//
	//	===== Autres Couches ====
	for (uint i=1; i < N; i++) {
		if (couche_type[i] == 2) {
			//	=== Couches Neurones ===
			mdl->couche_neurone_conn[i] = malloc(sizeof(uint*) * couche_y[i]);
			for (uint j=0; j < couche_y[i]; j++)
				mdl->couche_neurone_conn[i][j] = cpyuint(couche_neurone_conn[i][j], couche_n[i]);
			//
			mdl->poid_depart[i] = mdl->poids;
			mdl->conste_depart[i] = mdl->constantes;
			mdl->locd_depart[i] = mdl->locds;
			//
			mdl->poids += couche_y[i] * (2*couche_n[i] + 1);
			mdl->constantes += 0;
			mdl->locds += 1*couche_y[i];
		} else {
			// === Couches Filtriques ===
			mdl->couche_filtre_depart[i] = cpyuint(couche_filtre_depart[i], couche_y[i]);
			//
			mdl->poid_depart[i] = mdl->poids;
			mdl->conste_depart[i] = mdl->constantes;
			mdl->locd_depart[i] = mdl->locds;
			//
			mdl->poids += 0;
			mdl->constantes += couche_y[i]*couche_n[i];
			mdl->locds += couche_y[i]*(6+couche_n[i]+couche_n[i]-1);
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
	for (uint i=0; i < mdl->poids; i++) {
		mdl->poid[i] = (float)(rand()%10000)/9999.0 - 0.5;
	}
	for (uint i=0; i < mdl->constantes; i++) {
		mdl->constante[i] = (float)(rand()%10000)/9999.0;
	}
	//
	return mdl;
}

void liberer_mdl(Mdl_t * mdl) {
	uint N = mdl->N;
	//
	for (uint i=1; i < N; i++) {
		if (mdl->couche_type[i] == 2) {
			for (uint j=0; j < mdl->couche_y[i]; j++)
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

static char * noms[] = {"Filtre-Prix", "Filtres", "Neurone"};

void plume_mdl(Mdl_t * mdl) {
	printf("======== Mdl ==========\n");
	uint C = mdl->couches;
	for (uint i=0; i < C; i++) {
		printf("%i| %s :\n", i, noms[mdl->couche_type[i]]);
		for (uint j=0; j < mdl->couche_y[i]; j++) {
			if (mdl->couche_type[i] == 2) {
				printf("\ty[%i] : ", j);
				for (uint k=0; k < mdl->couche_n[i]; k++) {
					printf("%i, ", mdl->couche_neurone_conn[i][j][k]);
				}
				printf("\n");
			} else if (mdl->couche_type[i] == 1) {
				printf("\ty[%i] : depart = %i\n",
					j, mdl->couche_filtre_depart[i][j]);
			} else if (mdl->couche_type[i] == 0) {
				printf("\ty[%i] : ema=%i intervalle=%i\n",
					j, mdl->ema[j], mdl->intervalles[j]);
			}
		}
	}
};