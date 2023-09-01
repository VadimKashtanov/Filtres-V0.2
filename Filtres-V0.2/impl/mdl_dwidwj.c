#include "mdl.h"

uint inverse(uint n, float * tableau, float * inverse) {
	float coef;
	for (uint L=0; L < n; L++) {
		for (uint y=0; y < n; y++) {
			if (y == L) continue;
			//
			if (tableau[L*n + L] == 0) return 1;	//Erreur, matrice non inversible
			//
			coef = tableau[y*n + L]/tableau[L*n + L];
			for (uint k=0; k < n; k++) {
				tableau[y*n + k] -= tableau[L*n + k]*coef;
				inverse[y*n + k] -= inverse[L*n + k]*coef;
			}
		};
	};
	//
	for (uint L=0; L < n; L++) {
		coef = tableau[L*n + L];
		if (coef == 0) return 1;
			
		for (uint i=0; i < n; i++) {
			tableau[L*n + i] /= coef;
			inverse[L*n + i] /= coef;
		}
	}
	//
	return 0;
};