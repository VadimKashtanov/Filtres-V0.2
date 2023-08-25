#include "mdl.h"

#define P 0.05

static float estimer_alpha(Mdl_t * mdl, uint depart, uint N) {
	float somme = 0.0;
	//
	for (uint i=0; i < N; i++) {
		uint p0 = rand() % mdl->poids;
		uint p1 = rand() % mdl->poids;
		//
		float _f = f(mdl, depart);
		mdl->poid[p0] += 1e-5;
		float _fx = f(mdl, depart);
		mdl->poid[p1] += 1e-5;
		float _fxy = f(mdl, depart);
		mdl->poid[p0] -= 1e-5;
		float _fy = f(mdl, depart);
		mdl->poid[p1] -= 1e-5;
		//
		somme += (_fxy - _fx -_fy -_f)/1e-10;
	}
	//
	return fabs(1.0 / (somme != 0.0 ? somme : 1.0));
};

#define OPTI_TOUT_LES 1
#define ALPHA_TOUT_LES 10000

float score(Mdl_t * mdl) {
	uint nb_tests_second = (uint)roundf(1.0 + (float)mdl->poids*P);
	//
	float alpha = estimer_alpha(mdl, DEPART, nb_tests_second);
	//
	float gain;
	float score = 0;
	//
	float p1, p0;
	//
	for (uint i=DEPART; i < PRIXS-1; i++) {
		p1 = prixs[i+1];
		p0 = prixs[i];
		//
		gain = f(mdl, i) * (p1/p0 - 1) * LEVIER;
		score += gain;

		if (i % ALPHA_TOUT_LES == 0) {
			alpha = estimer_alpha(mdl, i, nb_tests_second);
		};

		if (i % OPTI_TOUT_LES == 0) {
			df(mdl, i, +gain);
			for (uint p=0; p < mdl->poids; p++)
				mdl->poid[p] -= alpha * mdl->d_poid[p];
		}
	}
	//
	return score / (PRIXS-DEPART);	//Score Moyen pour que je vois. Si besoin du vrai gain, utiliser printf()
};