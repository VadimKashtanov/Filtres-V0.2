#include "analyze_model.h"

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

	/*ptr(" === constantes ===\n");
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
	}*/

	/*ptr(" === grad ===\n");
	for (uint i=mdl->y_depart[1]; i < mdl->vars; i++)
		ptr("%i| %f\n", i, mdl->d_var[i]);

	*/ptr(" === var ===\n");
	for (uint i=0; i < mdl->vars; i++)
		ptr("%i| %f\n", i, mdl->var[i]);
};

void verifier_derivee(Mdl_t * mdl) {
	float _d_poid[mdl->poids];
	float _f = objectif_gain(mdl, DEPART);
	const float _1E5 = 1e-3;

	FOR(0, i, mdl->poids) {
		mdl->poid[i] += _1E5;
		_d_poid[i] = (objectif_gain(mdl, DEPART)-_f)/_1E5;
		mdl->poid[i] -= _1E5;
		mdl->d_poid[i] = 0;
	};

	d_objectif_gain(mdl, DEPART, objectif_gain(mdl, DEPART));
	
	FOR(0, i, mdl->poids) {
		float a = _d_poid[i];
		float b = mdl->d_poid[i];
		
		if (fabs(a+b)/2*.01 > fabs(a-b) || fabs(a-b) < 0.001) {
			printf("\033[42m%i|  %f  --  %f\033[0m\n", i, _d_poid[i], mdl->d_poid[i]);
		} else {
			printf("\033[41m%i|  %f  --  %f   [x%f]\033[0m\n",
				i, _d_poid[i], mdl->d_poid[i], _d_poid[i]/mdl->d_poid[i]);
		}
	};
};

//	===============================================================================
//	======================== Analyze comportement =================================
//	===============================================================================

void comportement(Mdl_t * mdl) {
	uint depart = DEPART + (rand() % (PRIXS-DEPART));

	float f_arr[30];

	for (uint i=0; i < 30; i++) {
		f_arr[i] = f(mdl, depart);
	};

	gnuplot(prixs + depart, 30, "Prixs");
	gnuplot(f_arr, 30, "Valeur de f (achat vente)");

	plume_mdl(mdl);
};