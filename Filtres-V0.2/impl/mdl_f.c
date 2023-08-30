#include "mdl.h"

/*
def filtre(x,f):
	n = len(x)
	#
	_min, _max = min(x), max(x)
	x = [(e-_min)/(_min-_max) for e in x]
	#
	s = sum((1+|xi-fi|)**.5 for i in 0..n) / n - 1
	d = sum((1+|xi1-xi  -  fi1+fi|)**2 for i in 0..(n-1)) / (n-1) -1 
	#
	return exp(-s*s -d*d)
*/

//	=======================================

static float filtre_n(float * locd, float * arr, float * filtre, uint n) {
	float _min=arr[0], _max=arr[0];
	float _x;
	float max_pos=0.0, min_pos=0.0;
	for (uint i=1; i < n; i++) {
		_x = arr[i];
		if (_x > _max) {
			_max = _x;
			max_pos = 1.0*i;
		}
		if (_x < _min) {
			_min = _x;
			min_pos = 1.0*i;
		}
	}

	locd[0] = _min;
	locd[1] = _max;
	locd[2] = max_pos;
	locd[3] = min_pos;

	//
	float x[n];

	//
	float _s = 0;
	const float d = _max-_min;
	float tmp, signe;
	for (uint i=0; i < n; i++) {
		x[i] = (arr[i]-_min)/d;
		//
		//
		tmp = x[i] - filtre[i];
		signe = (tmp >= 0 ? 1 : -1);
		tmp = sqrtf(1.0f+fabs(tmp));
		//
		_s += tmp;
		locd[4 + i] = 1/(2*tmp) * signe;
	}
	_s = _s/n - 1;

	//
	float _d = 0.0;
	for (int i=0; i < n-1; i++) {
		tmp = (x[i+1]-x[i]) - (filtre[i+1]-filtre[i]);
		signe = (tmp >= 0 ? 1 : -1); 
		tmp = 1.0f + fabs(tmp);
		//
		_d += powf(tmp,2);
		locd[4 + n + i] = 2*tmp*signe;
	}
	_d = _d/(n-1) - 1;

	locd[4 + n + n - 1] = _s;
	locd[4 + n + n -1 + 1] = _d;

	//
	return ___exp(-_s*_s - _d*_d);
};

static void df_filtre_n(
	//	Que du locd
	float * locd,
	//
	float dy, float * grad,
	float * arr, float * filtre, uint n)
{
	float _min=locd[0], _max=locd[1];
	uint _maxpos=locd[2], _minpos=locd[3];
	float *locd_xi_s=locd+4, *locd_i_d=locd+4+n;
	float _s=locd[4+n+(n-1)], _d=locd[4+n+(n-1)+1];
	//
	float dexp = ___exp(-_s*_s - _d*_d) * dy;	//dexp = exp
	float _d_s = -2*_s*dexp/n;
	float _d_d = -2*_d*dexp/(n-1);
	//
	float dx[n];
	//
	const float d = _max-_min;
	for (uint i=0; i < n; i++) {
		dx[i] = _d_s * locd_xi_s[i];
	}
	//
	for (uint i=0; i < n-1; i++) {
		dx[i+1] += _d_d*locd_i_d[i];
		dx[i] += -_d_d*locd_i_d[i];
	};
	//
	float _dmax=0, _dmin=0;
	for (uint i=0; i < n; i++) {
		grad[i] += dx[i]/d; 							
		grad[_maxpos] -= dx[i] * (arr[i]-_min)/(d*d);
		grad[_minpos] -= dx[i] / d;					
		grad[_minpos] += dx[i] * (arr[i]-_min)/(d*d);
	
	}
};
//	=======================================

static float neurone_n(float * locd, float * arr, float * poid, uint n) {
	float _somme = 0.0;
	float tmp;
	for (uint i=0; i < n; i++) {
		tmp = arr[i]*poid[i*2] + poid[i*2+1];
		locd[i] = tmp;
		_somme += ___tanh(tmp);
	}
	locd[n] = _somme;
	return ___tanh(_somme/n + poid[n*2]);	//ou gauss, a voire
};

static void d_neurone_n(
	float * locd,
	float dy, float * grad,
	float * arr, float * poid, float * d_poid, uint n)
{
	float _somme = locd[n];
	float _d_somme = ___d_tanh(_somme/n + poid[n*2]) * dy / n;
	d_poid[n*2] = _d_somme*n;
	//
	float tmp;
	for (uint i=0; i < n; i++) {
		tmp = ___d_tanh(locd[i]) * _d_somme;
		grad[i] = tmp * poid[i*2];
		d_poid[i*2] = tmp * arr[i];
		d_poid[i*2+1] = tmp;
	}
};

//	=======================================


float f(Mdl_t * mdl, uint depart) {
	uint C = mdl->N;
	uint * type = mdl->type;
	uint * y = mdl->y;
	uint * n = mdl->n;

	//	Filtres N
	float * x = allouer_flotants(n[0]);

	FOR(0, i, y[0]) {
		uint _ema = mdl->ema[i];
		uint interv = mdl->intervalles[i];
		//
		FOR(0, j, n[0]) {
			x[j] = ema[_ema][depart - j*interv];
		}

		mdl->var[i] = filtre_n(
			mdl->locd + mdl->locd_depart[0] + i*LOCDS_FLTR(n[0]), 
			x, mdl->constante + i*CONSTS_FLTR(n[0]), n[0]);

		mdl->d_var[i] = 0.0;
	};
	
	float * _x = allouer_flotants(mdl->max_n);
	FOR(1, i, C) {
		FOR(0, j, y[i]) {
			if (type[i] == 1) {
				mdl->var[mdl->y_depart[i] + j] = filtre_n(
					mdl->locd + mdl->locd_depart[i] + j*LOCDS_FLTR(mdl->n[i]),
					//
					mdl->var + mdl->y_depart[i-1] + mdl->fltr_depart[i][j],
					mdl->constante + mdl->conste_depart[i] + mdl->n[i]*j,
					mdl->n[i]
				);
			} else if (mdl->type[i] == 2) {
				FOR(0, k, n[i]) {
					_x[k] = mdl->var[mdl->y_depart[i-1] + mdl->neu_vers[i][j][k]];
				}
				//
				mdl->var[mdl->y_depart[i] + j] = neurone_n(
					mdl->locd + mdl->locd_depart[i] + j*LOCDS_NEU(mdl->n[i]),
					_x, mdl->poid + mdl->poid_depart[i] + j*(2*mdl->n[i]+1),
					mdl->n[i]
				);
			} else {
				ERR("Pas de couche %i", mdl->type[i]);
			}
			//
			mdl->d_var[mdl->y_depart[i] + j] = 0.0;	//pour pas faire 2 boucles
		}
	}
	free(_x);
	free(x);
	//
	return mdl->var[mdl->vars-1];
};

void df(Mdl_t * mdl, uint depart, float erreur) {
	mdl->d_var[mdl->vars-1] = erreur;

	uint C = mdl->N;
	uint * type = mdl->type;
	uint * y = mdl->y;
	uint * n = mdl->n;

	//
	//faire l'inverse
	//

	float * _x = allouer_flotants(mdl->max_n);
	float * _dx = allouer_flotants(mdl->max_n);

	RETRO_FOR(i, C, 1) {
		FOR(0, j, mdl->y[i]) {
			if (type[i] == 1) {
				FOR(0, k, n[i]) {
					_dx[k] = 0;
				}
				//
				df_filtre_n(
					mdl->locd + mdl->locd_depart[i] + j*LOCDS_FLTR(n[i]),
					//
					mdl->d_var[mdl->y_depart[i] + j],
					_dx,
					//
					mdl->var + mdl->y_depart[i-1] + mdl->fltr_depart[i][j],
					mdl->constante + mdl->conste_depart[i] + mdl->n[i]*j,
					mdl->n[i]
				);

				FOR(0, k, n[i]) {
					mdl->d_var[mdl->y_depart[i-1] + mdl->fltr_depart[i][j] + k] += _dx[k];
				}

			} else if (type[i] == 2) {
				FOR(0, k, n[i]) {
					_x[k] = mdl->var[mdl->y_depart[i-1] + mdl->neu_vers[i][j][k]];
					_dx[k] = 0.0;
				}
				
				d_neurone_n(
					mdl->locd + mdl->locd_depart[i] + j*LOCDS_NEU(n[i]),
					mdl->d_var[mdl->y_depart[i] + j],
					_dx,
					_x, mdl->poid + mdl->poid_depart[i] + j*POIDS_NEU(n[i]),
					mdl->d_poid + mdl->poid_depart[i] + j*POIDS_NEU(n[i]),
					mdl->n[i]
				);

				FOR(0, k, n[i]) {
					mdl->d_var[mdl->y_depart[i-1] + mdl->neu_vers[i][j][k]] += _dx[k];
				}
			}
		}
	}
	free(_x);
	free(_dx);
};

float objectif_gain(Mdl_t * mdl, uint depart) {
	float tmp = USDT * LEVIER * (prixs[depart+1]/prixs[depart]-1.0);
	return powf(f(mdl, depart)*tmp-fabs(tmp), 2)/2;
};

void d_objectif_gain(Mdl_t * mdl, uint depart, float obj_gain) {
	float tmp = USDT * LEVIER * (prixs[depart+1]/prixs[depart]-1.0);
	df(mdl, depart, -sqrtf(2*obj_gain)*tmp);
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
