#include "mdl.h"

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

	locd[0] = _max;
	locd[1] = _min;
	locd[2] = max_pos;
	locd[3] = min_pos;

	//
	float x[n];

	//
	float _s = 0;
	const float d = _max-_min;
	float tmp;
	for (int i=0; i < n; i++) {
		x[i] = (arr[i]-_min)/d;
		tmp = sqrtf(1.0f+fabs(x[i] - filtre[i]));
		_s += tmp;
		locd[4 + i] = 1/(2*tmp);
	}
	_s = _s/n - 1;

	//
	float _d = 0.0;
	float signe;
	for (int i=0; i < n-1; i++) {
		tmp = (x[i+1]-x[i]) - (filtre[i+1]-filtre[i]);
		signe = fabs(tmp)/tmp;
		tmp = 1.0f + fabs(tmp);
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
	float _min, float _max,
	uint _maxpos, uint _minpos,
	float * locd_xi_s,	//n
	float * locd_i_d,	//n-1
	float _s, float _d,
	//
	float dy, float * grad, float * arr, float * filtre, uint n)
{
	//
	//	==== return ___exp(-_s*_s - _d*_d);
	//
	float dexp = ___exp(-_s*_s - _d*_d) * dy;	//dexp = exp
	float _d_s = -2*_s*dexp/n;
	float _d_d = -2*_d*dexp/(n-1);
	//
	float dx[n];
	//
	const float d = _max-_min;
	for (uint i=0; i < n; i++) {
		dx[i] += _d_s * locd_xi_s[i];
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
		_dmax += (arr[i]-_min)/(d*d)*dx[i];
		_dmin += -dx[i]/d;
	}
	//
	grad[_maxpos] += _dmax;
	grad[_minpos] += _dmin;
};

//	=======================================

static float neurone_n(float * locd, float * arr, float * poid, uint n) {
	float _somme = 0.0;
	for (uint i=0; i < n; i++) _somme += ___tanh(arr[i]*poid[i*2] + poid[i*2+1]);
	locd[0] = _somme;
	return ___tanh(_somme/n + poid[n*2]);	//ou gauss, a voire
};

static void d_neurone_n(
	float _somme, float dy, float * grad,
	float * arr, float * poid, float * d_poid, uint n)
{
	float _d_somme = ___d_tanh(_somme/n + poid[n*2]) * dy / n;
	d_poid[n*2] = _d_somme*n;
	//
	float tmp;
	for (uint i=0; i < n; i++) {
		tmp = ___d_tanh(arr[i]*poid[i*2] + poid[i*2+1]) * _d_somme;
		grad[i] += tmp * poid[i*2];
		d_poid[i*2] = tmp * arr[i];
		d_poid[i*2+1] = tmp;
	}
};

//	=======================================


float f(Mdl_t * mdl, uint depart) {
	//	Filtres N
	uint n = mdl->couche_n[0];
	float x[n];
	for (uint i=0; i < mdl->couche_y[0]; i++) {
		uint _ema = mdl->ema[i];
		uint interv = mdl->intervalles[i];
		//
		for (uint j=0; j < n; j++) x[j] = ema[_ema][depart - j*interv];
		mdl->var[i] = filtre_n(
			mdl->locd + mdl->locd_depart[0] + i*(6+n*2-1), 
			x, mdl->constante + i*n, n);
	};
	//
	for (uint i=1; i < mdl->couches; i++) {
		for (uint j=0; j < mdl->couche_y[i]; j++) {
			float _x[mdl->couche_n[i]];
			//
			if (mdl->couche_type[i] == 1) {
				mdl->var[mdl->y_depart[i] + j] = filtre_n(
					mdl->locd + mdl->locd_depart[i] + j*(6+mdl->couche_n[i]*2-1),
					//
					mdl->var + mdl->y_depart[i-1] + mdl->couche_filtre_depart[i][j],
					mdl->constante + mdl->conste_depart[i] + mdl->couche_n[i]*j,
					mdl->couche_n[i]
				);
			} else if (mdl->couche_type[i] == 2) {
				for (uint k=0; k < mdl->couche_n[i]; k++)
					_x[k] = mdl->var[mdl->y_depart[i-1] + mdl->couche_neurone_conn[i][j][k]];
				mdl->var[mdl->y_depart[i] + j] = neurone_n(
					mdl->locd + mdl->locd_depart[i] + j,
					_x, mdl->poid + mdl->poid_depart[i] + j*(2*mdl->couche_n[i]+1),
					mdl->couche_n[i]
				);
			} else {
				ERR("Pas de couche %i", mdl->couche_type[i]);
			}
			//
			mdl->d_var[mdl->y_depart[i] + j] = 0.0;	//pour pas faire 2 boucles
		}
	}
	//
	return mdl->var[mdl->vars-1];
};

void df(Mdl_t * mdl, uint depart, float erreur) {
	mdl->d_var[mdl->vars-1] = erreur;
	//
	//faire l'inverse
	//
	for (int i=mdl->couches-1; i >= 1; i--) {
		
		for (uint j=0; j < mdl->couche_y[i]; j++) {
			
			float _x[mdl->couche_n[i]];
			float _dx[mdl->couche_n[i]];

			float * _locd;

			//
			if (mdl->couche_type[i] == 1) {
				_locd = mdl->locd + mdl->locd_depart[i] + j*(6+mdl->couche_n[i]*2-1);
				
				//mdl->var[mdl->y_depart[i] + j]
				df_filtre_n(
					_locd[0], _locd[1],
					(uint)_locd[2], (uint)_locd[3],
					_locd + 4, _locd + 4 + mdl->couche_n[i],
					_locd[4+mdl->couche_n[i]*2-1], _locd[4+mdl->couche_n[i]*2-1+1],
					//
					mdl->d_var[mdl->y_depart[i] + j],
					_dx,
					//
					mdl->var + mdl->y_depart[i-1] + mdl->couche_filtre_depart[i][j],
					mdl->constante + mdl->conste_depart[i] + mdl->couche_n[i]*j,
					mdl->couche_n[i]
				);

			} else if (mdl->couche_type[i] == 2) {
				for (uint k=0; k < mdl->couche_n[i]; k++) {
					_x[k] = mdl->var[mdl->y_depart[i-1] + mdl->couche_neurone_conn[i][j][k]];
					_dx[k] = mdl->d_var[mdl->y_depart[i-1] + mdl->couche_neurone_conn[i][j][k]];
				}
				d_neurone_n(
					mdl->locd[mdl->locd_depart[i] + j],
					mdl->d_var[mdl->y_depart[i] + j],
					_dx,
					_x, mdl->poid + mdl->poid_depart[i] + j*(2*mdl->couche_n[i]+1),
					mdl->d_poid + mdl->poid_depart[i] + j*(2*mdl->couche_n[i]+1),
					mdl->couche_n[i]
				);
				for (uint k=0; k < mdl->couche_n[i]; k++) {
					mdl->d_var[mdl->y_depart[i-1] + mdl->couche_neurone_conn[i][j][k]] = _dx[k];
				}
			}
		}
	}
};