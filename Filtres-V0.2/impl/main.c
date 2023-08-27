#include "mdl.h"

Mdl_t * generer(uint C, uint * type, uint * y, uint * n) {
	assert(C > 0);
	assert(type[0] == 0);
	assert(y[C-1] == 1);
	//
	uint *** couche_neurone_conn = malloc(sizeof(uint**) * C);
	uint ** couche_filtre_depart = malloc(sizeof(uint*) * C);
	//
	couche_filtre_depart[0] = 0;
	couche_neurone_conn[0] = 0;
	//
	for (uint i=1; i < C; i++) {
		if (type[i] == 1) {
			couche_neurone_conn[i] = 0;
			couche_filtre_depart[i] = malloc(sizeof(uint) * y[i]);
			for (uint j=0; j < y[i]; j++) {
				couche_filtre_depart[i][j] = rand() % (y[i-1]-n[i]+1);
			}
		} else if (type[i] == 2) {
			couche_filtre_depart[i] = 0;
			couche_neurone_conn[i] = malloc(sizeof(uint*) * y[i]);
			for (uint j=0; j < y[i]; j++) {
				couche_neurone_conn[i][j] = malloc(sizeof(uint) * n[i]);
				for (uint k=0; k < n[i]; k++) {
					couche_neurone_conn[i][j][k] = rand() % y[i-1];
				}
			};
		} else {
			ERR("type[i] == %i", type[i]);
		}
	};
	//
	uint _intervalles[y[0]];
	uint _ema[y[0]];
	for (uint i=0; i < y[0]; i++) {
		_intervalles[i] = intervalles[rand() % INTERVALLES];
		_ema[i] = rand() % NB_DIFF_EMA;
	}
	//
	Mdl_t * mdl = cree_mdl(
		C,
		type, n, y,
		couche_neurone_conn,
		couche_filtre_depart,
		_intervalles,
		_ema);
	//
	for (uint i=1; i < C; i++) {
		if (type[i] == 1) {
			free(couche_filtre_depart[i]);
		} else if (type[i] == 2) {
			for (uint j=0; j < y[i]; j++)
				free(couche_neurone_conn[i][j]);
			free(couche_neurone_conn[i]);
		} else {
			ERR("type[i] == %i", type[i]);
		}
	};
	free(couche_filtre_depart);
	free(couche_neurone_conn);
	//
	return mdl;
};

#define M 1000

/*void selection() {
	float meilleur_score;
	//
	for (uint i=0; i < M; i++) {
		//
		mdl = meilleur_mdl
		//
		if (score > meilleur_score) {
			meilleur_score = score;
			ecrire_model(mdl, "model.bin");
		}
	}
};*/

void verifier_derivee(Mdl_t * mdl) {
	float _d_poid[mdl->poids];
	float _f = objectif_gain(mdl, DEPART);
	const float _1E5 = 1e-3;
	for (uint i=0; i < mdl->poids; i++) {
		mdl->poid[i] += _1E5;
		_d_poid[i] = (objectif_gain(mdl, DEPART)-_f)/_1E5;
		mdl->poid[i] -= _1E5;
		mdl->d_poid[i] = 0;
	};
	//

	//	!!!!!!!!!!!!!!!
	//	Il faut grad(gain) et pas grad(f).
	//  !!!!!!!!!!!!!!!

	d_objectif_gain(mdl, DEPART, objectif_gain(mdl, DEPART));
	//
	for (uint i=0; i < mdl->poids; i++) {
		float a = _d_poid[i];
		float b = mdl->d_poid[i];
		//
		if (fabs(a+b)/2*.01 > fabs(a-b) || fabs(a-b) < 0.0001) {
			printf("\033[42m%i|  %f  --  %f\033[0m\n", i, _d_poid[i], mdl->d_poid[i]);
		} else {
			printf("\033[41m%i|  %f  --  %f\033[0m\n", i, _d_poid[i], mdl->d_poid[i]);
		}
	};
};

#define NEU 2
#define FLTR 1

int main() {
	srand(0);
	charger_les_prixs();
	//
	uint couche_type[] = { 0, NEU, FLTR, NEU};
	uint couche_y[] =    { 4,   4,    3,   1};
	uint couche_n[] =    { 6,   2,    3,   3};
	//
	Mdl_t * mdl = generer(4, couche_type, couche_y, couche_n);
	//
	verifier_derivee(mdl);
	plume_mdl(mdl);
	//
	//printf("%f\n", f(mdl, DEPART));
	//raise(SIGINT);
	//
	liberer_mdl(mdl);
};