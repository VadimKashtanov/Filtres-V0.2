#include "mdl.h"

#include "analyze_model.h"
#include "score.h"

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
		_intervalles[i] = intervalles[0];//rand() % INTERVALLES];
		_ema[i] = 0;//rand() % NB_DIFF_EMA;
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

typedef struct {
	uint c, y, n;
} A_t;

Mdl_t * gen(A_t * l, uint taille) {
	uint couche_type[taille];
	uint y[taille];
	uint n[taille];
	for (uint i=0; i < taille; i++) {
		couche_type[i] = l[i].c;
		y[i] = l[i].y;
		n[i] = l[i].n;
	}
	return generer(taille, couche_type, y, n);
};

void mise_a_jour_seconde(Mdl_t * mdl, uint depart) {
	uint n = mdl->poids;
	float * tableau = allouer_flotants(n*n);
	float * _inverse = allouer_flotants(n*n);
	for (uint i=0; i < n; i++) for (uint j=0; j < n; j++) _inverse[i*n+j] = 1.0*(uint)(i==j);
	//
	float passe[n];
	
	d_objectif_gain(mdl, depart, objectif_gain(mdl, depart));
	memcpy(passe, mdl->d_poid, sizeof(float) * n);
	//
	for (uint i=0; i < n; i++) {
		mdl->poid[i] += 1e-3;
		d_objectif_gain(mdl, depart, objectif_gain(mdl, depart));
		for (uint j=0; j < n; j++)
			tableau[i*n + j] = (mdl->d_poid[j] - passe[j])/1e-3 + (rnd()*0.0001);
		mdl->poid[i] -= 1e-3;
	}

	if (inverse(n, tableau, _inverse) != 0) return;

	float somme;
	for (uint i=0; i < n; i++) {
		somme = 0;
		for (uint k=0; k < n; k++) {
			somme += 0.1*_inverse[i*n + k] * mdl->d_poid[k]/n;
		}
		mdl->poid[i] += somme;
	}	
};

static float * filtre_alpha_mdl(Mdl_t * mdl, float * les_alpha) {
	float * ret = malloc(sizeof(float) * mdl->poids);
	//
	FOR(0, c, mdl->C) {
		if (mdl->type[c] == 2) {
			FOR(0, j, POIDS_NEU(mdl->n[c])*mdl->y[c]) {
				ret[mdl->poid_depart[c] + j] = les_alpha[c];
			}
		}
	}
	//
	return ret;
};

int main() {
	srand(0);
	charger_les_prixs();

	MODE_OBJECTIF = 0;

#define N 3

	A_t pile[N] = {
		{.c=0,    .y=3, .n=6},
		//{.c=FLTR, .y=32, .n=4},
		//
		//{.c=NEU,  .y=16, .n=6},
		//{.c=NEU,  .y=8, .n=4},
		{.c=NEU,  .y=2, .n=2},
		{.c=NEU,  .y=1, .n=2}
	};
	float les_alpha[] = {
		0,
		//0,
		//0.1,
		//0.00001,
		0.2,
		0.1
	};
	
	Mdl_t * mdl = gen(pile, N);
	float * filtre_alpha = filtre_alpha_mdl(mdl, les_alpha);

	verifier_derivee(mdl);

	zero_dpoid(mdl);
	derivee_et_seconde(mdl, DEPART + (rand() % PRIXS-DEPART-1));
	//zero_dpoid(mdl);
	//derivee_et_seconde(mdl, DEPART + (rand() % PRIXS-DEPART-1));

	//comportement(mdl);
	//comportement(mdl);
	//comportement(mdl);

	//for (uint i=0; i < 100; i++) {
		//FOR(0, p, mdl->poids) mdl->poid[p] = 2*rnd()-1;
		//score(mdl, les_alpha);
	//}
	srand(0);
	comportement(mdl);
	printf("##################################\n");

	score(mdl, filtre_alpha);
	score(mdl, les_alpha);
	score(mdl, les_alpha);

	zero_dpoid(mdl);
	derivee_et_seconde(mdl, DEPART + (rand() % PRIXS-DEPART-1));

	srand(0);
	comportement(mdl);
	//derivee_et_seconde(mdl, DEPART + (rand() % PRIXS-DEPART-1));

	/*//printf("%f\n", estimer_alpha(mdl, DEPART, 5));

	srand(0);
	comportement(mdl);
	//
	//for (uint i=0; i < 100; i++)
	//	mise_a_jour_seconde(mdl, DEPART + (rand() % PRIXS-DEPART-1));
	//
	//getchar();
	score(mdl);

	srand(0);
	comportement(mdl);

	uint n = mdl->poids;

	printf("%f\n", dp2(mdl, DEPART, 0, 0));
	printf("%f\n", dp2(mdl, DEPART, 10, 10));
	printf("%f\n", dp2(mdl, DEPART, n-1, n-1));

	printf("%f\n", dp(mdl, DEPART, 0));
	printf("%f\n", dp(mdl, DEPART, 10));
	printf("%f\n", dp(mdl, DEPART, n-1));*/
	
	//
	liberer_mdl(mdl);
};