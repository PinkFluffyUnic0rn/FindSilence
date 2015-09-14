#include <stdint.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

typedef struct _complexd
{
	double real;
	double imag;
} complexd;

int log_2( uint32_t a );

complexd *fft( complexd *data, uint samples_count, int inv );
