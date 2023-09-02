#include "analyze_model.h"

static char * noms[] = {"Filtre-Prix", "Filtres", "Neurone"};

static void plume_poids(Mdl_t * mdl) {
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
};

static void plume_constes(Mdl_t * mdl) {
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
};

static void plume_var(Mdl_t * mdl) {
	ptr(" === vars ===\n");
	FOR(0, i, mdl->C) {
		ptr("#%i %s\n", i, noms[mdl->type[i]]);
		FOR(0, j, mdl->y[i]) {
			ptr("%i| %f \n", mdl->y_depart[i] + j, mdl->var[mdl->y_depart[i] + j]);
		}
	}
};

static void plume_grad(Mdl_t * mdl) {
	ptr(" === grad ===\n");
	FOR(0, i, mdl->C) {
		ptr("#%i %s\n", i, noms[mdl->type[i]]);
		FOR(0, j, mdl->y[i]) {
			ptr("%i| %f \n", mdl->y_depart[i] + j, mdl->d_var[mdl->y_depart[i] + j]);
		}
	}
};

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

	plume_poids(mdl);
	plume_constes(mdl);
	plume_var(mdl);
	plume_grad(mdl);
};

void verifier_derivee(Mdl_t * mdl) {
	uint depart = DEPART + (rand()%(PRIXS-DEPART-1));
	//
	float _d_poid[mdl->poids];
	float _f = objectif_gain(mdl, depart);
	const float _1E5 = 1e-4;

	FOR(0, i, mdl->poids) {
		mdl->poid[i] += _1E5;
		_d_poid[i] = (objectif_gain(mdl, depart)-_f)/_1E5;
		mdl->poid[i] -= _1E5;
		mdl->d_poid[i] = 0;
	};

	d_objectif_gain(mdl, depart, objectif_gain(mdl, depart));
	
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
#define T 50

	uint depart = DEPART + (rand() % (PRIXS-DEPART-T-1));

	float f_arr[T];

	float var_par_t[T*mdl->vars];

	for (uint i=0; i < T; i++) {
		f_arr[i] = f(mdl, depart + i);
		memcpy(var_par_t + i*mdl->vars, mdl->var, sizeof(float)*mdl->vars);
	};

	for (uint i=0; i < mdl->vars; i++) {
		//printf("%i| ", i);
		for (uint t=0; t < T; t++) {
		//	if (var_par_t[t*mdl->vars + i] >= 0) printf(" ");
		//	printf("%f | ", var_par_t[t*mdl->vars + i]);
		}
		//printf("\n");
	}

	//gnuplot(ema[1] + depart - 6*6, 6, "ema1");
	//gnuplot(ema[1] + depart + T - 6*6, 6, "ema2");

	//gnuplot(prixs + depart, T, "Prixs");
	gnuplot(f_arr, T, "Valeur de f (achat vente)");

	//plume_mdl(mdl);
	//plume_poids(mdl);
	//plume_constes(mdl);
};

//================================================================================

#include "score.h"

void derivee_et_seconde(Mdl_t * mdl, uint depart) {
	printf(" === derivee et derivee second sur : depart=%i === \n", depart);
	//
	d_objectif_gain(mdl, depart, objectif_gain(mdl, depart));
	//
	for (uint c=0; c < mdl->couches; c++) {
		printf("v -- %s -- v\n", noms[mdl->type[c]]);
		if (mdl->type[c] == 2) {
			for (uint j=0; j < POIDS_NEU(mdl->n[c])*mdl->y[c]; j++) {
				uint i = mdl->poid_depart[c] + j;
				printf("%3.i| dw = %s%f\n",//  dwidwi=%s%f\n",
					i,
					(mdl->d_poid[i] >= 0 ? " " : ""), mdl->d_poid[i]//,
					//(dp2(mdl, DEPART, i, i) >= 0 ? " " : ""), dp2(mdl, DEPART, i, i)
				);
			}
		}
	};
};