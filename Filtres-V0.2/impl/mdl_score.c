#include "score.h"

#define P 0.05

float dp2(Mdl_t * mdl, uint depart, uint p0, uint p1) {
	static const float _1E5 = 1e-3;
	//
	float _f = objectif_gain(mdl, depart);
	mdl->poid[p0] += _1E5;
	float _fx = objectif_gain(mdl, depart);
	mdl->poid[p1] += _1E5;
	float _fxy = objectif_gain(mdl, depart);
	mdl->poid[p0] -= _1E5;
	float _fy = objectif_gain(mdl, depart);
	mdl->poid[p1] -= _1E5;
	//printf("%f %f %f %f\n", _fxy, _fx, _fy, _f);
	return (_fxy - _fx -_fy -_f)/(_1E5*_1E5);
};

float dp(Mdl_t * mdl, uint depart, uint p) {
	static const float _1E5 = 1e-3;
	//
	float _f = objectif_gain(mdl, depart);
	mdl->poid[p] += _1E5;
	float _fx = objectif_gain(mdl, depart);
	mdl->poid[p] -= _1E5;
	//printf("%f %f %f %f\n", _fxy, _fx, _fy, _f);
	return (_fx-_f)/(_1E5);
};

float estimer_alpha(Mdl_t * mdl, uint depart, uint N) {
	float somme = 0.0;
	//
	for (uint i=0; i < N; i++) {
		uint p0 = rand() % mdl->poids;
		uint p1 = rand() % mdl->poids;
		//
		somme += dp2(mdl, depart, p0, p1);
	}
	//printf("%f\n", somme);
	//printf("somme = %f\n", somme);
	//
	float alpha = fabs(somme/N);
	//
	if (somme != 0) alpha = 1.0 / fabs(somme/N); 
	if (somme == 0) alpha = 1.0;
	//
	return alpha;
};

#define OPTI_TOUT_LES 1000
#define ALPHA_TOUT_LES 100000000000

float score(Mdl_t * mdl, float * les_alpha) {
	uint nb_tests_second = (uint)roundf(1.0 + (float)mdl->poids*P);
	//printf("%i\n", nb_tests_second);
	//
	float alpha = estimer_alpha(mdl, DEPART, nb_tests_second);
	//
	float _gain, _score;
	float gain_total=0;
	float score = 0;
	//
	float p1, p0;
	memset(mdl->d_poid, 0, sizeof(float) * mdl->poids);


	float * suivie = allouer_flotants((PRIXS-DEPART)/100);

UNE_COURBE(dp0);

	//
	for (uint i=DEPART; i < PRIXS-1; i++) {
		p1 = prixs[i+1];
		p0 = prixs[i];
		//
		_score = objectif_gain(mdl, i);
		score += _score;
		//
		_gain = mdl->var[mdl->vars - 1] * (p1/p0 - 1) * LEVIER;
		gain_total += _gain;
		if ((i-DEPART)%100==0) suivie[(i-DEPART)/100] = gain_total;
		//

		if (i % ALPHA_TOUT_LES == 0) {
			//alpha = estimer_alpha(mdl, i, nb_tests_second);
			//
			//printf("%i/%i\n", i, PRIXS);
		};

		d_objectif_gain(mdl, i, _score);
		if (i % OPTI_TOUT_LES == 0) {
			for (uint p=0; p < mdl->poids; p++) {
				//printf("%f\n", mdl->d_poid[p]);
				//if (mdl->d_poid[p] != 0) mdl->poid[p] += _score/mdl->d_poid[p];
				mdl->poid[p] -= /*alpha*/ les_alpha[p]*mdl->d_poid[p] / (float)OPTI_TOUT_LES;// / (1+p*10);
				SUIVIE_COURBE(dp0, mdl->d_poid[0]);
				mdl->d_poid[p] = 0;
			}
			//printf("=============================\n");
		}
	}
	PLUMER_LA_COURBE(dp0);
	//
	//gnuplot(suivie, (PRIXS-DEPART)/100, "suivie des gains");
	printf("Gain total = %f\n", gain_total);
	//
	return score / (PRIXS-DEPART);	//Score Moyen pour que je vois. Si besoin du vrai gain, utiliser printf()
};