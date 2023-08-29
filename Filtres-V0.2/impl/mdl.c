#include "mdl.h"

Mdl_t * cree_mdl(
	uint C,
	uint * type, uint * n, uint * y,
	uint *** neu_vers,	//[i] = 0 si i-eme est non neuronal
	uint ** fltr_depart,	//[i] = 0 si ieme est non filtrique
	// 0 -eme couche
	uint * _intervalles,
	uint * _ema)
{
	assert(C > 0);
	assert(type[0] == 0);
	assert(y[C-1] == 1);
	assert(n[0] <= MAX_N_PREMIERE_COUCHE);
	//
	Mdl_t * mdl = malloc(sizeof(Mdl_t));
	//
	uint N = C;
	mdl->N = C;
	mdl->C = C;
	mdl->couches = C;
	mdl->type = cpyuint(type, C);
	mdl->n = cpyuint(n, C);
	mdl->y = cpyuint(y, C);
	//
	mdl->max_n = u_max(n, C);
	//
	mdl->intervalles = cpyuint(_intervalles, y[0]);
	mdl->ema = cpyuint(_ema, y[0]);
	//
	for (uint i=0; i < y[0]; i++) {
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
	mdl->neu_vers = calloc(C, sizeof(uint**));
	mdl->fltr_depart = calloc(C, sizeof(uint*));
	//
	//	==== Couche 0 =====
	mdl->poid_depart[0] = 0;
	mdl->conste_depart[0] = 0;
	mdl->locd_depart[0] = 0;
	mdl->y_depart[0] = 0;
	//
	mdl->vars = y[0];
	mdl->poids = 0;
	mdl->constantes = y[0] * CONSTS_FLTR(n[0]);
	mdl->locds = y[0] * LOCDS_FLTR(n[0]);
	//
	//	===== Autres Couches ====
	FOR(1, i, N) {
		if (mdl->type[i] == 2) {
			//	=== Couches Neurones ===
			mdl->neu_vers[i] = malloc(sizeof(uint*) * y[i]);
			
			FOR(0, j, y[i]) {
				mdl->neu_vers[i][j] = cpyuint(neu_vers[i][j], n[i]);
			}
			//
			mdl->poid_depart[i] = mdl->poids;
			mdl->conste_depart[i] = mdl->constantes;
			mdl->locd_depart[i] = mdl->locds;
			//
			mdl->poids += y[i] * POIDS_NEU(n[i]);
			mdl->constantes += 0;
			mdl->locds += y[i] * LOCDS_NEU(n[i]);
		} else {
			// === Couches Filtriques ===
			mdl->fltr_depart[i] = cpyuint(fltr_depart[i], y[i]);
			//
			mdl->poid_depart[i] = mdl->poids;
			mdl->conste_depart[i] = mdl->constantes;
			mdl->locd_depart[i] = mdl->locds;
			//
			mdl->poids += 0;
			mdl->constantes += y[i] * CONSTS_FLTR(n[i]);
			mdl->locds += y[i] * LOCDS_FLTR(n[i]);
		}
		//
		mdl->y_depart[i] = mdl->vars;
		mdl->vars += y[i];
	}
	//
	mdl->var = allouer_flotants(mdl->vars);
	mdl->d_var = allouer_flotants(mdl->vars);
	mdl->d_poid = allouer_flotants(mdl->poids);
	mdl->constante = allouer_flotants(mdl->constantes);
	mdl->poid = allouer_flotants(mdl->poids);
	mdl->locd = allouer_flotants(mdl->locds);
	//
	FOR(0, i, mdl->poids) {
		mdl->poid[i] = 2*((float)(rand()%10000)/9999.0 - 0.5);
	}

	FOR(0, i, mdl->constantes) {
		mdl->constante[i] = (float)(rand()%10000)/9999.0;
	}
	//
	return mdl;
}

void liberer_mdl(Mdl_t * mdl) {
	uint N = mdl->N;
	//
	for (uint i=1; i < N; i++) {
		if (mdl->type[i] == 2) {
			FOR(0, j, mdl->y[i]) {
				free(mdl->neu_vers[i][j]);
			}
			free(mdl->neu_vers[i]);
		} else {
			free(mdl->fltr_depart[i]);

		}
	}
	free(mdl->neu_vers);
	free(mdl->fltr_depart);
	//
	free(mdl->n);
	free(mdl->y);
	free(mdl->type);
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

void ecrire_mdl(Mdl_t * mdl, char * fichier) {
	FILE * fp = fopen(fichier, "wb");
	fwrite(&mdl->N, sizeof(uint), 1, fp);
	uint C = mdl->N;
	//
	fwrite(mdl->type, sizeof(uint), C, fp);
	fwrite(mdl->y, sizeof(uint), C, fp);
	fwrite(mdl->n, sizeof(uint), C, fp);
	//
	fwrite(mdl->intervalles, sizeof(uint), mdl->y[0], fp);
	fwrite(mdl->ema, sizeof(uint), mdl->y[0], fp);
	//
	FOR(1, i, C) {
		if (mdl->type[i] == 2) {
			FOR(0, j, mdl->y[i]) {
				fwrite(mdl->neu_vers[i][j], sizeof(uint), mdl->n[i], fp);
			}
		} else {
			fwrite(mdl->fltr_depart[i], sizeof(uint), mdl->n[i], fp);
		}
	}
	//
	fwrite(mdl->constante, sizeof(float), mdl->constantes, fp);
	fwrite(mdl->poid, sizeof(float), mdl->poids, fp);
	//
	fclose(fp);
};

Mdl_t * lire_mdl(char * fichier) {
	FILE * fp = fopen(fichier, "rb");
	//
	uint C;
	fread(&C, sizeof(uint), 1, fp);
	//
	uint y[C], n[C], type[C];
	//
	fread(type, sizeof(uint), C, fp);
	fread(y, sizeof(uint), C, fp);
	fread(n, sizeof(uint), C, fp);
	//
	uint intervalles[C], ema[C];
	fread(intervalles, sizeof(uint), y[0], fp);
	fread(ema, sizeof(uint), y[0], fp);
	//
	uint ** neu_vers[C];
	uint * fltr_depart[C];
	for (uint i=1; i < C; i++) {
		if (type[i] == 2) {
			neu_vers[i] = malloc(sizeof(uint*) * y[i]);
			
			FOR(0, j, y[i]) {
				neu_vers[i][j] = malloc(sizeof(uint) * n[i]);
				fread(neu_vers[i][j], sizeof(uint), n[i], fp);
			}

		} else {
			fltr_depart[i] = malloc(sizeof(uint) * n[i]);
			fread(fltr_depart[i], sizeof(uint), n[i], fp);
		}
	}
	//
	Mdl_t * mdl = cree_mdl(
		C,
		type, n, y,
		neu_vers,	//[i] = 0 si i-eme est non neuronal
		fltr_depart,	//[i] = 0 si ieme est non filtrique
		intervalles,
		ema);
	//
	fread(mdl->constante, sizeof(float), mdl->constantes, fp);
	fread(mdl->poid, sizeof(float), mdl->poids, fp);
	//
	fclose(fp);
	//
	for (uint i=1; i < C; i++) {
		if (type[i] == 2) {
			FOR(0, j, mdl->y[i]) {
				free(neu_vers[i][j]);
			}
			free(neu_vers[i]);
		} else {
			free(fltr_depart[i]);
		}
	}
	//
	return mdl;
};

static char * noms[] = {"Filtre-Prix", "Filtres", "Neurone"};

void plume_mdl(Mdl_t * mdl) {
	printf("======== Mdl ==========\n");
	uint C = mdl->couches;
	for (uint i=0; i < C; i++) {
		printf("%i| %s n=%i:\n", i, noms[mdl->type[i]], mdl->n[i]);
		for (uint j=0; j < mdl->y[i]; j++) {
			if (mdl->type[i] == 2) {
				printf("\ty[%i] : ", j);
				for (uint k=0; k < mdl->n[i]; k++) {
					printf("%i, ", mdl->neu_vers[i][j][k]);
				}
				printf("\n");
			} else if (mdl->type[i] == 1) {
				printf("\ty[%i] : depart = %i\n",
					j, mdl->fltr_depart[i][j]);
			} else if (mdl->type[i] == 0) {
				printf("\ty[%i] : ema=%i intervalle=%i\n",
					j, mdl->ema[j], mdl->intervalles[j]);
			}
		}
	};

	ptr(" === poids ===\n");
	FOR(0, i, mdl->C) {
		ptr("%s #%i\n", noms[mdl->type[i]], i);
		if (mdl->type[i] == 2) {
			uint poids = POIDS_NEU(mdl->n[i]);
			FOR(0, j, mdl->y[i]) {
				ptr(" .y%i:\n", j);
				FOR(0, k, poids) {
					ptr("%i| %f\n",
					mdl->poid_depart[i] + j*poids + k,
					mdl->poid[mdl->poid_depart[i] + j*poids + k]);
				}
			}
		}
	}

	ptr(" === constantes ===\n");
	FOR(0, i, mdl->C) {
		ptr("%s #%i\n", noms[mdl->type[i]], i);
		if (mdl->type[i] < 2) {
			uint constantes = CONSTS_FLTR(mdl->n[i]);
			FOR(0, j, mdl->y[i]) {
				ptr(" .y%i:\n", j);
				FOR(0, k, constantes) {
					ptr("%i| %f\n",
					mdl->conste_depart[i] + j*constantes + k,
					mdl->constante[mdl->conste_depart[i] + j*constantes + k]);
				}
			}
		}
	}

	ptr(" === grad ===\n");
	for (uint i=mdl->y_depart[1]; i < mdl->vars; i++)
		ptr("%i| %f\n", i, mdl->d_var[i]);

	ptr(" === var ===\n");
	for (uint i=0; i < mdl->vars; i++)
		ptr("%i| %f\n", i, mdl->var[i]);
};