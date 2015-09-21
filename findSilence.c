#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

#include "fft.h"

#define MIN_FREQUENCY (200 * 2)
#define MAX_FREQUENCY (4000 * 2)

#define SECS_IN_PART 3.0
#define MEANS_IN_PART 10

struct audio_file
{
	void *raw_data;
	uint samples_count;
	uint frequency;
	uint depth;
	uint channels;
};

#ifdef DEBUG
#include <time.h>

#include <SDL/SDL.h>
#include "findSilence_test.h"
#endif

void print_event( uint c, uint i, char ev, uint frequency )
{
	uint t = (double) i * 1000.0 / (double) frequency;
	printf( "%u %u %c\n", c, t, ev );
}

void print_output( const uint *silence, uint len, uint c,
	const struct audio_file *a_file )
{
	char type = 'q';
	char prev_type;
	uint i;
	
	print_event( c, 0, 'b', a_file->frequency );
	
	for ( i = 0; i < len; ++i )
	{
		prev_type = type;
		type = silence[i] ? 's' : 'S';

		if ( prev_type != type )
		{
			uint cur_sample = (double) i
				* (double) (a_file->samples_count)
				/ (double) len;

			print_event( c, cur_sample, type, a_file->frequency );
		}
	}

	print_event( c, a_file->samples_count, 'e', a_file->frequency );
}


double derivite_sqr( const double *sample_data, uint beg, uint end )
{
	double s = 0;
	uint i;

	for ( i = beg; (i + 1) <= end; ++i )
	{
		double a = (sample_data[i+1]-sample_data[i]);

		s += a * a;
	}
	
	return s / pow( (double)(end - beg), 2.0 );
}

void signal_derivitive_sqr( double **sample_data, uint len, uint d )
{
	uint newlen = len / d + ((len % d) ? 1 : 0);
	double *d_sig = (double *) malloc( sizeof(double) * newlen );
	uint i, j;

	j = 0;
	for ( i = 0; i < len; i += d )
	{
		if ( (i + d) > len )
			d = len - i;
		
		d_sig[j++] = derivite_sqr( *sample_data, i, i + d );	
	}

	free( *sample_data );

	*sample_data = d_sig;
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
		sample_data[i] /= abs(max-min) / 1000.0;
	
	return sample_data;
}

double average( const double *sample_data, uint b, uint e )
{
	double s = 0;
	uint i;

	for ( i = b; i < e; ++i )
		s += sample_data[i];

	return s / (e - b);
}

void smooth_signal( double **in_data, uint samples_count, uint d )
{
	double *out_data = (double *) malloc( sizeof(double)
		* samples_count / d );
	uint i;

	for ( i = 0; (i+d) < samples_count; i += d )
		out_data[i/d] = average( *in_data, i, i + d );

	free( *in_data );

	*in_data = out_data;
}

double *get_channel( const struct audio_file *a_file, uint c )
{
	double *s_data
		= (double *) malloc( a_file->samples_count * sizeof(double) );
	uint i;
	
	for ( i = 0; i < a_file->samples_count; ++i )
		switch( a_file->depth )
		{
		case 1:
			s_data[i] = ((int8_t *) a_file->raw_data)[i*a_file->channels + c];
			break;
		case 2:
			s_data[i] = ((int16_t *) a_file->raw_data)[i*a_file->channels + c];
			break;
		}
	
	return s_data; 
}

void load_wav( const char *fname, struct audio_file *a_file )
{
	FILE *file;
	char buf[255];
	uint pcm = 0;
	
	memset( a_file, 0, sizeof(struct audio_file) );

	file = fopen( fname, "r" );
	if ( file == NULL )
	{
		fprintf( stderr, "Cannot load file\n" ); 
		exit( 0x02 );
	}

	fseek( file, 0, SEEK_END );
	if ( ftell(file) < 44 )
	{
		fprintf( stderr, "File format not supported\n" ); 
		exit( 0x12 );
	}
	fseek( file, 0, SEEK_SET );


	fread( buf, sizeof(char), 4, file );	
	if ( strncmp(buf, "RIFF", 4) )
	{
		fprintf( stderr, "File format not supported\n" ); 
		exit( 0x22 );
	}

	fseek( file, 4, SEEK_CUR );
	
	fread( buf, sizeof(char), 4, file );
	if ( strncmp(buf, "WAVE", 4) )
	{
		fprintf( stderr, "File format not supported\n" ); 
		exit( 0x32 );
	}

	fread( buf, sizeof(char), 4, file );
	if ( strncmp(buf, "fmt", 3) )
	{
		fprintf( stderr, "File format not supported\n" ); 
		exit( 0x42 );
	}

	fseek( file, 4, SEEK_CUR );

	fread( &pcm, sizeof(char), 2, file );
	if ( pcm != 1 )
	{
		fprintf( stderr, "File format not supported\n" ); 
		exit( 0x52 );
	}

	fread( &(a_file->channels), sizeof(char), 2, file );
	
	fread( &(a_file->frequency), sizeof(char), 4, file );
	
	fseek( file, 4, SEEK_CUR );
	
	fread( &(a_file->depth), sizeof(char), 2, file );

	fseek( file, 2, SEEK_CUR );
	
	fread( buf, sizeof(char), 4, file );
	if ( strncmp(buf, "data", 4) )
	{
		fprintf( stderr, "File format not supported\n" ); 
		exit( 0x62 );
	}

	fread( &(a_file->samples_count), sizeof(char), 4, file );

	a_file->raw_data = malloc( a_file->samples_count );
	
	fread( a_file->raw_data, sizeof(char), a_file->samples_count, file );

	a_file->samples_count /= a_file->depth;
	a_file->depth /= a_file->channels;

	fclose( file );
}

double metric( double a, double b )
{
	return pow( a - b, 2.0 );
}

void assign_centroids( const double *sample_data, uint *sample_centroids,
	uint len, const double *means, uint k )
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

int recompute_centroids( const double *sample_data,
	const uint *sample_centroids, uint len, double *means, uint k )
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

void init_centroids( const uint *hist, double *means,
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
	for ( i = min; i < max / 2; ++i )
		s0 += hist[i];

	s1 = 0;
	for ( i = max / 2; i < max; ++i )
		s1 += hist[i];

	p0 = (double) s0 / (double) (s0+s1);
	p1 = (double) s1 / (double) (s0+s1);

	k0 = floor(p0 * (double) k + 0.5);
	k1 = floor(p1 * (double) k + 0.5);

	if ( (double)((max + min) / 2 - min) / (double) k0 > 3.0
		&& (double)(max - (max + min) / 2) / (double) k1 > 3.0 )
	{
		init_centroids( hist, means, k_min, k_min + k0, min,
			(max + min) / 2 );
		init_centroids( hist, means, k_min + k0, k_max,
			(max + min) / 2, max );
	}
	else
	{
		double d = (double)(max - min) / (double)(k_max - k_min);

		means[k_min] = min;
		for ( i = k_min + 1; i < k_max; ++i )
			means[i] = means[i - 1] + d;

		return;
	}
}

uint *k_means( const double *sample_data, uint len, uint k, double **means )
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

uint *find_silence( const double *sample_data, uint len, uint part_len )
{
	uint k = MEANS_IN_PART;
	uint *silence = malloc( sizeof(uint) * len );
	uint i, j;

	for ( i = 0; i < len; i += part_len )
	{
		double *means;
		uint *clusters;
		uint id_min = 0;

		if ( (i + part_len) >= len )
		{
			k = k * (len - i) / part_len;
			part_len = len - i;
		}
		
		clusters = k_means( sample_data + i, part_len, k, &means );

		for ( j = 0; j < k; ++j )
			id_min = (means[j] < means[id_min]) ? j : id_min;

		for ( j = 0; j < part_len; ++j )
			silence[i + j] = (clusters[j] == 0);
	
		free( clusters );
	}

	return silence;
}

void remove_short_ranges( uint *silence, uint len,
	uint min_len_si, uint val )
{
	uint *tmp_pr, *tmp;

	for ( tmp = silence; tmp < (silence + len); ++tmp )
	{
		tmp_pr = tmp;
		uint c;
	
		c = 0;
		while ( *tmp == val && tmp < (silence + len) )
		{
			++tmp;
			++c;
		}
			
		tmp = (tmp < (silence + len)) ? tmp : (silence + len);

		if ( c < min_len_si )
			for ( ; tmp_pr != tmp; ++tmp_pr )
				*tmp_pr = !val;	
	}
}

double *calculate_frequency_sd( const double *sample_data,
	uint samples_count, uint d, uint win_size )
{
	uint i, j;
		
	double *sdf = (double *) malloc( sizeof(complexd)
		* samples_count / d ); 
		
	double *weights = (double *) malloc( sizeof(double) * win_size );
	for ( i = 0; i < win_size; ++i )
		weights[i] = 0.5 * ( 1.0 - cos(2.0 * M_PI *
				(double) i / (double) (win_size - 1)) );

	complexd *arg = (complexd *) malloc( sizeof(complexd)
		* win_size );
	complexd *c  = (complexd *) malloc( sizeof(complexd)
		* win_size );
	
	for ( j = 0; j < samples_count; j += d )
	{
		for ( i = 0; i < win_size; ++i )
		{
			int idx = j + i - ( win_size / 2 );
			if ( (idx < samples_count) && (idx >= 0) )
				arg[i].real = (double) sample_data[idx]
					* weights[i];
			else
				arg[i].real = 0;
	
		}
	
		fft( arg, c, win_size, 1 );
		
		double mean = 0.0, sd = 0.0;
	
		for ( i = 0; i < win_size; ++i )
			mean += sqrt( pow(c[i].real, 2.0)
				+ pow(c[i].imag, 2.0) );
		
		mean /= win_size;
	
		for ( i = 0; i < win_size; ++i )
			sd += pow( mean - sqrt( pow(c[i].real, 2.0)
				+ pow(c[i].imag, 2.0)), 2.0 );
		
		sd /= win_size;
		sd = sqrt( sd );
		
		sdf[j / d] = sd;
		
	}

	free(weights);
	free(arg);
	free(c);

	return sdf;
}

int main( int argc, const char **argv )
{
	struct audio_file a_file;
	double border;
	double min_len_si_secs;
	double min_len_so_secs;
	double *sfm;
	int i, j;

	load_wav( argv[1], &a_file );
	border = atof( argv[2] );
	min_len_si_secs = atof( argv[3] );
	min_len_so_secs = atof( argv[4] );

	
	for ( i = 0; i < a_file.channels; ++i )
	{

// for testing //
#ifdef DEBUG
		struct timespec t0, t1;

		clock_gettime( CLOCK_REALTIME, &t0 );		
#endif

		double *sample_data = get_channel( &a_file, i );
		uint len = a_file.samples_count;
		uint d;
		uint *silence = malloc( sizeof(uint) * len );
	

		to_borders( sample_data, len );

		d = a_file.frequency / MAX_FREQUENCY;
		smooth_signal( &sample_data, len, d );
		len /= d;	

		d = 1;
		signal_derivitive_sqr( &sample_data, len, d );
		len /= d;

		sfm = calculate_frequency_sd( sample_data, len,
			MAX_FREQUENCY / MIN_FREQUENCY, 512 );

		len = len * MIN_FREQUENCY / MAX_FREQUENCY; 

		free( sample_data );

		sample_data = sfm;
	
		if ( border > 0.0 )
		{
			d = MIN_FREQUENCY / 10;
			signal_derivitive_sqr( &sample_data, len, d );
			len /= d;

			silence = malloc( sizeof(uint) * len );
			for ( j = 0; j < len; ++j )
				silence[j] = (sample_data[j] < border) ? 1 : 0;	
		}
		else
		{
			uint part_len = SECS_IN_PART * a_file.frequency
				* len / a_file.samples_count;
			
			silence = find_silence( sample_data, len, part_len );
		}
	
		uint min_len_si = min_len_si_secs * (double) (a_file.frequency)
			* (double) len / (double) (a_file.samples_count);
		uint min_len_so = min_len_so_secs * (double) (a_file.frequency)
			* (double) len / (double) (a_file.samples_count);

		remove_short_ranges( silence, len, min_len_si, 1 );
		remove_short_ranges( silence, len, min_len_so, 0 );

		print_output( silence, len, i, &a_file );

// for testing //
#ifdef DEBUG
		clock_gettime( CLOCK_REALTIME, &t1 );

		long long int dt = t1.tv_sec * 1000000000LL + t1.tv_nsec
			- t0.tv_sec * 1000000000LL + t0.tv_nsec;
		
		printf( "time elapsed: %f sec\n", (double) dt / 1e9 );


		if ( SDL_Init(SDL_INIT_AUDIO) < 0 )
		{
			fprintf( stderr, "Cannot initilize SDL: %s\n",
				SDL_GetError() );
			exit( 1 );
		}

		draw_signal( sample_data, len, silence );
		draw_histogram( sample_data, len );

		silence_to_zero( &a_file, silence, len );

		play_audio( &a_file, a_file.samples_count*2 );
		
		fgetc(stdin);

		SDL_CloseAudio();
#endif
	}

	free( a_file.raw_data );

	return 0;
}
