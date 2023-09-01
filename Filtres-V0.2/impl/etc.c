#include "etc.h"

float rnd() {
	return (float)(rand()%100)/99.0;
};

inline float ___exp(register float x)  // cubic spline approximation
{
    union { float f; int i; } reinterpreter;

    reinterpreter.i = (int)(12102203.0f*x) + 127*(1 << 23);
    int m = (reinterpreter.i >> 7) & 0xFFFF;  // copy mantissa
    // empirical values for small maximum relative error (8.34e-5):
    reinterpreter.i +=
         ((((((((1277*m) >> 14) + 14825)*m) >> 14) - 79749)*m) >> 11) - 626;
    return reinterpreter.f;
}

inline float ___gauss(register float x) {return ___exp(-x*x);};
inline float ___d_gauss(register float x) {return -2*x*___gauss(x);};

inline float ___logistique(register float x) {return 1.0/(1.0+___exp(-x));};//2*___tanh(x)+0.5;};    //  2*(tanh(x))+0.5
inline float ___d_logistique(register float x) {return ___logistique(x)*(1 - ___logistique(x));};

inline float ___tanh(register float x) {return tanhf(x);};//x/(0.5 + fabs(x));};		//  x/( 0.5 + fabs(x) )
inline float ___d_tanh(register float x) {return 1 - powf(___tanh(x), 2);};

void gnuplot(float * arr, uint len, char * titre) {
	char buff[200];
	//
	FILE * fp = fopen("gnuplot_dat.dat", "w");
	//
	for (uint i=0; i < len; i++) {
		snprintf(buff, 100, "%i ", i);
		fputs(buff, fp);
		//
		snprintf(buff, 100, "%f\n", arr[i]);
		fputs(buff, fp);
	}
	fclose(fp);
	//
	snprintf(
		buff,
		200,
		"gnuplot -p -e \"set title \'%s\'; plot 'gnuplot_dat.dat' w lp\"",
		titre);
	//
	assert(!system(buff));
	//
	assert(!system("rm gnuplot_dat.dat"));
};

uint* cpyuint(uint * arr, uint len) {
	uint * ret = malloc(sizeof(uint) * len);
	memcpy(ret, arr, sizeof(uint) * len);
	return ret;
}

float* allouer_flotants(uint nb) {
	return malloc(sizeof(float) * nb);
}

uint u_max(uint * x, uint len) {
	uint _max=x[0];
	FOR(1,i,len) {
		if (x[i] > _max)
			_max = x[i];
	};
	return _max;
}