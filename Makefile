all:
	gcc findSilence.c -o findsilence \
	`pkg-config --cflags --libs gtk+-3.0 allegro_acodec-5.0  allegro-5.0` -lm
