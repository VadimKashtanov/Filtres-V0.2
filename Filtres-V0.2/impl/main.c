#include "mdl.h"

#define M 1000

void selection() {
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
};

int main() {
	charger_les_prixs();
	//
	//
	uint couche_type[] = {0,1,2,2,2,2};
	uint couche_y[] = {64, 16, 32, 8, 4, 1};
	uint couche_n[] = {6, 4, 8, 8, 4, 4};
	//
	uint couche_neurone_conn[6] = {};
	uint couche_filtre_depart[6] = {};
	//
	uint _intervalles[64] = {1};
	uint _ema[64] = {0};
	//
	Mdl_t * mdl = cree_mdl(
		6,
		couche_type, couche_n, couche_y,
		couche_neurone_conn,
		couche_filtre_depart,
	);
	//
	float _d_poid[];
	//
};