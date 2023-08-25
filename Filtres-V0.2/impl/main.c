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
				couche_filtre_depart[i][j] = rand() % (y[i-1]-n[i]);
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
	float _f = f(mdl, DEPART);
	for (uint i=0; i < mdl->poids; i++) {
		mdl->poid[i] += 1e-5;
		_d_poid[i] = (f(mdl, DEPART)-_f)/1e-5;
		mdl->poid[i] -= 1e-5;
		mdl->d_poid[i] = 0;
	};
	//

	//	!!!!!!!!!!!!!!!
	//	Il faut grad(gain) et pas grad(f).
	//  !!!!!!!!!!!!!!!

	df(mdl, DEPART, f(mdl, DEPART));
	//
	for (uint i=0; i < mdl->poids; i++) {
		if (fabs(_d_poid[i]-mdl->d_poid[i]) < 0.0001) {
			printf("\033[42m%i|  %f  --  %f\033[0m\n", i, _d_poid[i], mdl->d_poid[i]);
		} else {
			printf("\033[41m%i|  %f  --  %f\033[0m\n", i, _d_poid[i], mdl->d_poid[i]);
		}
	};
};

int main() {
	srand(0);
	charger_les_prixs();
	//
	uint couche_type[] = {0,1,2,2,2,2};
	uint couche_y[] = {64, 16, 32, 8, 4, 1};
	uint couche_n[] = {6, 4, 8, 8, 4, 4};
	//
	Mdl_t * mdl = generer(6, couche_type, couche_y, couche_n);
	//
	verifier_derivee(mdl);
	//
	//printf("%f\n", f(mdl, DEPART));
	//raise(SIGINT);
	//
	liberer_mdl(mdl);
};