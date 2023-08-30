#include "mdl.h"

Mdl_t * generer(uint C, uint * type, uint * y, uint * n) {
	assert(C > 0);
	assert(type[0] == 0);
	assert(y[C-1] == 1);
	//
	uint *** neu_vers = malloc(sizeof(uint**) * C);
	uint ** fltr_depart = malloc(sizeof(uint*) * C);
	//
	fltr_depart[0] = 0;
	neu_vers[0] = 0;
	//
	FOR(1, i, C) {
		if (type[i] == 1) {
			neu_vers[i] = 0;
			fltr_depart[i] = malloc(sizeof(uint) * y[i]);
			
			FOR(0, j, y[i]) {
				fltr_depart[i][j] = rand() % (y[i-1]-n[i]+1);
			}
		} else if (type[i] == 2) {
			fltr_depart[i] = 0;
			neu_vers[i] = malloc(sizeof(uint*) * y[i]);

			FOR(0, j, y[i]) {
				neu_vers[i][j] = malloc(sizeof(uint) * n[i]);
				//
				FOR(0, k, n[i]) {
					neu_vers[i][j][k] = rand() % y[i-1];
				}
			};
		} else {
			ERR("type[i] == %i", type[i]);
		}
	};
	//
	uint _intervalles[y[0]];
	uint _ema[y[0]];
	
	FOR(0, i, y[0]) {
		_intervalles[i] = intervalles[rand() % INTERVALLES];
		_ema[i] = rand() % NB_DIFF_EMA;
	}
	
	Mdl_t * mdl = cree_mdl(
		C,
		type, n, y,
		neu_vers,
		fltr_depart,
		_intervalles,
		_ema);
	//
	FOR(1, i, C) {
		if (type[i] == 1) {
			free(fltr_depart[i]);
		} else if (type[i] == 2) {
			FOR(0, j, y[i]) {
				free(neu_vers[i][j]);
			}
			free(neu_vers[i]);
		} else {
			ERR("type[i] == %i", type[i]);
		}
	};
	free(fltr_depart);
	free(neu_vers);
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

#define NEU 2
#define FLTR 1

int main() {
	srand(0);
	charger_les_prixs();

	//
	uint couche_type[] = {  0,   NEU, FLTR, NEU};
	uint y[] =    		 {  3,     3,    2,   1};
	uint n[] =    		 {  4,     3,    3,   2};

	//
	Mdl_t * mdl = generer(4, couche_type, y, n);
	ecrire_mdl(mdl, "mdl0");
	liberer_mdl(mdl);

	mdl = lire_mdl("mdl0");
	
	//
	liberer_mdl(mdl);
};