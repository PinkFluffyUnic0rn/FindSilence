#include "fft.h"

int log_2( uint32_t a )
{
	int c = 0;

	while ( a >>= 1 ) 
		++c;

	return c;
}

complexd *fft( const complexd *data, complexd *out, uint samples_count, int inv )
{
	int bits_count = log_2( samples_count );
	uint32_t i;
	int j;
	int k;

	for ( i = 0; i < samples_count; ++i )
	{
		int b = 0;
		int a = i;
			
		for ( j = 0; j < bits_count; ++j )
		{
			b <<= 1;
			b |= (a & 1);
			a >>= 1;
		}
	
		out[i] = data[b];
	}
	
	int lp;
	int lp2;
	for ( k = 1; k <= bits_count; ++k )
	{
		lp = pow( 2, k );
		lp2 = lp / 2;
		double darg = -inv * M_PI / lp2;
		double arg = 0;

		for ( j = 0; j < lp2; ++j )
		{
			double c = cos(arg);
			double s = sin(arg);

			arg += darg;

			for ( i = j; i < samples_count; i += lp )
			{
				int iw = i + lp2;
				double wr = out[iw].real * c - out[iw].imag * s;
				double wi = out[iw].real * s + out[iw].imag * c;
				
				out[iw].real = out[i].real - wr;
				out[iw].imag = out[i].imag - wi;
				out[i].real = out[i].real + wr;
				out[i].imag = out[i].imag + wi;
			}
		}
	}

	if ( inv == 1 )
		for ( i = 0; i < samples_count; ++i )
		{
			out[i].real = out[i].real / samples_count;
			out[i].imag = out[i].imag / samples_count;
		}
		
	return out;
}


void fft2d( complexd *data, int w, int h, int inv )
{
/*
	complexd *tmph = (complexd *) malloc( sizeof(complexd) * h );
	complexd *out;
	int i, j;

	for ( i = 0; i < h; ++i )
	{	
		out = fft( data + i*w, w, inv );
		
		memmove( data + i*w, out, sizeof(complexd) * w );
		free( out );
	}


	for ( i = 0; i < w; ++i )
	{
		for ( j = 0; j < h; ++j )
			tmph[j] = data[j * w + i];
		
		out = fft( tmph, h, inv );

		for ( j = 0; j < h; ++j )
			data[j * w + i] = out[j];

		free( out );
	}

	free( tmph );
*/
}


