#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

#include <allegro5/allegro.h>
#include <allegro5/allegro_audio.h>
#include <allegro5/allegro_acodec.h> 

#include "findSilence_test.h"

#define MIN_FREQUENCY 200
#define MAX_FREQUENCY 4000

#define SECS_IN_PART 5
#define MEANS_IN_PART 10

double border;
double min_length_secs;

char *fname;
ALLEGRO_SAMPLE *audio;

void *file_data;
uint samples_count = 0;
uint audio_frequency = 0;
ALLEGRO_AUDIO_DEPTH audio_depth = 0;
uint audio_channels = 0;

void print_event( uint c, uint i, char ev )
{
	printf( "%u %u %c\n", c, i*1000/audio_frequency, ev );
}

void print_output( uint *silence, uint len, uint c )
{
	char type = 'q';
	char prev_type;
	uint i;

	print_event( c, 0, 'b' );
	
	for ( i = 0; i < len; ++i )
	{
		prev_type = type;
		type = silence[i] ? 's' : 'S';

		if ( prev_type != type )
			print_event( c, i*samples_count/len, type );
	}

	print_event( c, samples_count, 'e' );
}

double derivite_sqr( double *sample_data, uint beg, uint end )
{
	double s = 0;
	double avr = 0.0f;
	uint i;

	for ( i = beg; (i + 1) < end; ++i )
	{
		double a = (sample_data[i+1]-sample_data[i]);

		s += a * a;
	}
	
	return s / pow( (double)(end - beg), 2.0 );
}
	
double *signal_derivitive_sqr( double *sample_data, uint len, uint d )
{
	uint newlen = len / d + ((len % d) ? 1 : 0);
	double *d_sig = (double *) malloc( sizeof(double) * newlen );
	uint i, j;

	j = 0;
	for ( i = 0; i < len; i += d )
	{
		if ( (i + d) >= len )
			d = len - i;
		
		d_sig[j++] = derivite_sqr( sample_data, i, i + d );
	}

	free( sample_data );

	return d_sig;
}

double *to_borders( double *sample_data, uint len )
{
	double max = sample_data[0];
	double min = max;
	uint i;

	for ( i = 0; i < len; ++i )
	{
		max = (sample_data[i] > max) ? sample_data[i] : max;
		min = (sample_data[i] < min) ? sample_data[i] : min;
	}

	for ( i = 0; i < len; ++i )
		sample_data[i] /= abs(max-min) / 1000;
	
	return sample_data;
}

double average( double *sample_data, uint b, uint e )
{
	double s = 0;
	uint i;

	for ( i = b; i < e; ++i )
		s += sample_data[i];

	return s / (e - b);
}

double *smooth_signal( double *in_data, uint samples_count, uint d )
{
	double *out_data = (double *) malloc( sizeof(double)
		* samples_count / d );
	uint i;

	for ( i = 0; (i+d) < samples_count; i += d )
		out_data[i/d] = average( in_data, i, i + d );

	free( in_data );

	return out_data;
}

double *get_channel( int16_t *sample_data, uint c )
{
	double *s_data = (double *) malloc( samples_count * sizeof(double) );
	uint i;

	for ( i = 0; i < samples_count; ++i )
		switch( audio_depth )
		{
		case ALLEGRO_AUDIO_DEPTH_INT8:
			s_data[i] = ((int8_t *) file_data)[i*audio_channels + c];
			break;
		case ALLEGRO_AUDIO_DEPTH_INT16:
			s_data[i] = ((int16_t *) file_data)[i*audio_channels + c];
			break;
		}

	return s_data; 
}

int load_file()
{
	audio = al_load_sample( fname );

	if ( !audio )
	{
		fprintf( stderr, "File format not supported\n" );
		exit( 1 );	
	}

	file_data = al_get_sample_data( audio );

	samples_count = al_get_sample_length( audio );
	audio_depth = al_get_sample_depth( audio );
	audio_frequency = al_get_sample_frequency( audio );
	audio_channels = al_get_channel_count( al_get_sample_channels(audio) );

	if ( audio_depth != ALLEGRO_AUDIO_DEPTH_INT8
		&& audio_depth != ALLEGRO_AUDIO_DEPTH_INT16 )
	{
		fprintf( stderr, "PCM not supported\n" );
		exit( 2 );
	}

	return 1;
}

void read_args( int argc, char **argv )
{
	fname = argv[1];
	border = atof( argv[2] );
	min_length_secs = atof( argv[3] );
}

void init_allegro()
{
	al_init();
	al_install_audio();
	al_init_acodec_addon();

	al_reserve_samples( 1 );
}

double metric( double a, double b )
{
	return pow( a - b, 2.0 );
}

void assign_centroids( double *sample_data, uint *sample_centroids, uint len,
	double *means, uint k )
{
	int i, j;

	for ( i = 0; i < len; ++i )
	{
		double min_dist = metric( sample_data[i], means[0] );
		int min_id = 0;

		for ( j = 1; j < k; ++j )
		{
			double dist = metric( sample_data[i], means[j] );
		
			if ( dist < min_dist )
			{
				min_dist = dist;
				min_id = j;
			}

		}

		sample_centroids[i] = min_id;
	}
}

int recompute_centroids( double *sample_data, uint *sample_centroids, uint len,
	double *means, uint k )
{
	double *s = malloc( sizeof(double) * k );
	uint *n = malloc( sizeof(double) * k );
	int i;

	for ( i = 0; i < k; ++i )
	{
		s[i] = 0.0;
		n[i] = 0;
	}

	for ( i = 0; i < len; ++i )
	{
		s[sample_centroids[i]] += sample_data[i];
		++n[sample_centroids[i]];
	}

	int stop = 1;
	for ( i = 0; i < k; ++i )
	{
		if ( n[i] != 0 )
		{
			s[i] /= n[i];	

			if ( abs( s[i] - means[i] ) > 0.001f )
				stop = 0;

			means[i] = s[i];
		}
		else
			means[i] = INFINITY;
	}

	free(s);
	free(n);

	return stop;
}

void init_centroids( uint *hist, double *means,
	uint k_min, uint k_max, uint min, uint max )
{
	uint s0, s1;
	uint i;
	double p0, p1;
	uint k0, k1;
	uint k = k_max - k_min;

	if ( k == 0 )
		return;

	s0 = 0;
	for ( i = min; i < max/2; ++i )
		s0 += hist[i];

	s1 = 0;
	for ( i = max/2; i < max; ++i )
		s1 += hist[i];

	p0 = (double)s0/(double)(s0+s1);
	p1 = (double)s1/(double)(s0+s1);

	k0 = floor(p0*(double)k + 0.5);
	k1 = floor(p1*(double)k + 0.5);

//	if ( (double)((max+min)/2 - min) / (double)(k0) > 25.0
//		&& (double)(max - (max+min)/2) / (double)(k1) > 25.0 )

	if ( (double) s0 / (double)(k0) > 50.0
		&& (double)((max+min)/2 - min) / (double)(k0) > 50.0
		&& (double) s1 / (double)(k1) > 50.0 
		&& (double)(max - (max+min)/2) / (double)(k1) > 50.0 )
	{
		init_centroids( hist, means, k_min, k_min + k0, min, (max+min)/2 );
		init_centroids( hist, means, k_min + k0, k_max, (max+min)/2, max );
	}
	else
	{
		double d = (double)(max - min) / (double)(k_max - k_min);

		means[k_min] = min;
		for ( i = k_min + 1; i < k_max; ++i )
			means[i] = means[i-1] + d;

		return;
	}

}

uint *k_means( double *sample_data, uint len, uint k, double **means )
{
	*means = malloc( sizeof(double) * k );
	uint *sample_centroids = malloc( sizeof(double) * len );
	uint *hist;
	uint i;

	double max = sample_data[0];
	double min = max;

	memset( sample_centroids, 0, sizeof(uint) * len );

	if ( k == 0 )
		return sample_centroids;

	for ( i = 0; i < len; ++i )
	{
		max = (sample_data[i] > max) ? sample_data[i] : max;
		min = (sample_data[i] < min) ? sample_data[i] : min;
	}

	hist = malloc( sizeof(uint) * ((uint)max + 1) );

	for ( i = 0; i < max; ++i )
		hist[i] = 0;

	for ( i = 0; i < len; ++i )
		++hist[(uint) abs(sample_data[i])];

	init_centroids( hist, *means,
		0, k, min, max + 1 );

	do
	{
		assign_centroids( sample_data, sample_centroids, len,
			*means, k );
	} while ( !recompute_centroids( sample_data, sample_centroids, len,
		*means, k ) );

	return sample_centroids;
}

uint *find_silence( double *sample_data, uint len )
{
	uint k = MEANS_IN_PART;
	uint d = SECS_IN_PART * audio_frequency * len / samples_count;
	uint *silence = malloc( sizeof(uint) * len );
	uint i, j;

	for ( i = 0; i < len; i += d )
	{
		double *means;
		uint *clusters;
		uint id_min = 0;

		if ( (i + d) >= len )
		{
			k = k * (len - i) / d;
			d = len - i;
		}
		
		clusters = k_means( sample_data + i, d, k, &means );
/*
for ( j = 0; j < k; ++j )
	printf( "* %f\n", means[j] );
printf( "\n" );
*/
		for ( j = 0; j < k; ++j )
			id_min = (means[j] < means[id_min]) ? j : id_min;

		for ( j = 0; j < d; ++j )
			silence[i + j] = (clusters[j] == 0);
	
		free( clusters );
	}

	return silence;
}

void remove_short_ranges( uint *silence, uint len )
{
	uint *tmp_pr, *tmp;
	uint min_length_samples = min_length_secs
		* (double)(audio_frequency * len / samples_count);
	
	for ( tmp = silence; ; ++tmp )
	{
		tmp_pr = tmp;
		uint c = 0;

		while ( *tmp != 0 && tmp < (silence + len) )
		{
			++tmp;
			++c;
		}
			
		if ( tmp >= (silence + len) )
			break;

		if ( c < min_length_samples )
			for ( ; tmp_pr != tmp; ++tmp_pr )
				*tmp_pr = 0;	
	}
}

int main( int argc, char **argv )
{
	int i, j;
	
	init_allegro();
	
	read_args( argc, argv );

	load_file();

	for ( i = 0; i < audio_channels; ++i )
	{
		double *sample_data = get_channel( file_data, i );
		uint len = samples_count;
		uint d;
		uint *silence;

		d = audio_frequency / MAX_FREQUENCY;
		sample_data = smooth_signal( sample_data, len, d );
		len /= d;	

		to_borders( sample_data, len );

		d = audio_frequency / d / MIN_FREQUENCY;
		sample_data = signal_derivitive_sqr( sample_data, len, d );
		len /= d;
		
		if ( border > 0.0 )
		{
			d = MIN_FREQUENCY/10;
			sample_data = signal_derivitive_sqr( sample_data, len, d );
			len /= d;

			silence = malloc( sizeof(uint) * len );
			for ( j = 0; j < len; ++j )
				silence[j] = (sample_data[j] < border) ? 1 : 0;	
		}
		else
			silence = find_silence( sample_data, len );

		remove_short_ranges( silence, len );
		
		print_output( silence, len, i );

// for testing //
		draw_signal( sample_data, len, silence );
		draw_histogram( sample_data, len );
	
		int16_t *raw_data = (int16_t *) file_data;
		
		for ( i = 0; i < samples_count; ++i )
			if ( silence[i/(samples_count/len)] )
				raw_data[i*audio_channels] = 0;
	
		printf( "%s\n", "EOF to stop playing" );
		al_play_sample( audio, 1.0f, 0.0f, 1.0f,
			ALLEGRO_PLAYMODE_ONCE, NULL );
		while ( fgetc(stdin) != EOF ) {}
	}

	al_destroy_sample( audio );

	return 0;
}
