#ifndef FIND_SILENCE_TEST_H
#define FIND_SILENCE_TEST_H

#include <cairo.h>
#include <gtk/gtk.h>

typedef unsigned int uint;

struct rgb
{
	unsigned char b, g, r, a;
};

struct audio_play
{
	Uint8 *audio_pos;
	Uint32 audio_len;
};

#define WIDTH 1024
#define HEIGHT 768

void draw_signal( double *sample_data, uint len, uint *silence )
{
	cairo_surface_t *csur;
	struct rgb *image_data;
	struct rgb blue = { 255, 0, 0, 255 }, red = { 0, 0, 255, 255 };
	double max, min;
	int i, j;	

	csur = cairo_image_surface_create( CAIRO_FORMAT_RGB24, WIDTH, HEIGHT );
	image_data = ((struct rgb *) cairo_image_surface_get_data( csur ));

	max = min = sample_data[0];
	for ( i = 0; i < len; ++i )
	{
		max = (sample_data[i] > max) ? sample_data[i] : max;
		min = (sample_data[i] < min) ? sample_data[i] : min;
	}

	for ( i = 0; i < HEIGHT; ++i )
		for ( j = 0; j < WIDTH; ++j )
		{
			image_data[i * WIDTH + j].a = 255;
			image_data[i * WIDTH + j].r = 0;
			image_data[i * WIDTH + j].g = 0;
			image_data[i * WIDTH + j].b = 0;
		}

	for ( i = 0; i < len; ++i )
	{
		int val = HEIGHT - 1 - sample_data[i] * HEIGHT / (max - min);
		int y = val;
		int x = (double)((WIDTH-1) * i) / (double)(len);
		
		y = (y > HEIGHT) ? (HEIGHT - 1) : ((y < 0) ? 0 : y);
	
		for ( j = y; j < HEIGHT; ++j )
			image_data[j*WIDTH+x] = silence[i] ? blue : red;	
	}

	cairo_surface_write_to_png( csur, "Signal" );
}

void draw_histogram( double *sample_data, uint len )
{
	cairo_surface_t *csur;
	struct rgb *image_data;
	uint *hist;
	uint i, k;
	double max, min;
	uint hist_size;

	csur = cairo_image_surface_create( CAIRO_FORMAT_RGB24, WIDTH, HEIGHT );
	image_data = ((struct rgb *) cairo_image_surface_get_data( csur ));

	max = min = sample_data[0];
	for ( i = 0; i < len; ++i )
	{
		max = (sample_data[i] > max) ? sample_data[i] : max;
		min = (sample_data[i] < min) ? sample_data[i] : min;
	}
	hist_size = max + 1;

	hist = calloc( hist_size, sizeof(uint) );
	for ( i = 0; i < len; ++i )
		++hist[(uint) abs(sample_data[i])];

	for ( i = 0; i < HEIGHT; ++i )
		for ( k = 0; k < WIDTH; ++k )
		{
			image_data[i*WIDTH+k].a = 255;
			image_data[i*WIDTH+k].r = 0;
			image_data[i*WIDTH+k].g = 0;
			image_data[i*WIDTH+k].b = 0;
		}

	for ( i = 0; i <= max; ++i )
	{
		uint y = hist[i];
		uint x = i * WIDTH / hist_size;
		y = (y > HEIGHT) ? (HEIGHT - 1) : ((y < 0) ? 0 : y);
		
		y = HEIGHT - y;
		for ( k = HEIGHT - 1; k > y; --k )
		{
			image_data[k*WIDTH+x].a = 255;
			image_data[k*WIDTH+x].r = 0;
			image_data[k*WIDTH+x].g = 255;
			image_data[k*WIDTH+x].b = 0;
		}
	}
	
	cairo_surface_write_to_png( csur, "Histogram" );	
}

void silence_to_zero( struct audio_file *a_file, uint *silence, uint len )
{
	int i;

	for ( i = 0; i < a_file->samples_count; ++i )
		if ( silence[i/(a_file->samples_count/len)] )
			switch ( a_file->depth )
			{
			case 1:
				((int8_t *) (a_file->raw_data))[i*a_file->channels] = 0;
				break;
			case 2:
				((int16_t *) (a_file->raw_data))[i*a_file->channels] = 0;
				break;
			}
}

void my_audio_callback( void *userdata, Uint8 *stream, int len )
{
	struct audio_play *play = (struct audio_play *) userdata;	

	len = ( len > play->audio_len ? play->audio_len : len );
	memcpy( stream, play->audio_pos, len );
	
	play->audio_pos += len;
	play->audio_len -= len;
}

void play_audio( struct audio_file *a_file, uint32_t len )
{
	struct audio_play *play = malloc( sizeof(struct audio_play) );	
	SDL_AudioSpec audio_spec;

	audio_spec.freq = a_file->frequency;

	switch ( a_file->depth )
	{
	case 1:
		audio_spec.format = AUDIO_S8;
		break;
	case 2:
		audio_spec.format = AUDIO_S16;
		break;
	}
	audio_spec.channels = a_file->channels;
	audio_spec.samples = 4096;
	audio_spec.callback = my_audio_callback;
	audio_spec.userdata = play;
	
	play->audio_pos = a_file->raw_data;
	play->audio_len = len;
	
	if ( SDL_OpenAudio(&audio_spec, NULL) < 0 )
	{
		fprintf( stderr, "Couldn't open audio: %s\n", SDL_GetError() );
		exit(-1);
	}
	
	SDL_PauseAudio(0);
}

#endif
