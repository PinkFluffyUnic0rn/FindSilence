all: release

release:
	gcc findSilence.c fft.c -o findsilence -Wall -lm

debug:
	gcc findSilence.c fft.c -o findsilence -Wall -g -pg -D DEBUG \
	`pkg-config --cflags --libs gtk+-3.0` -lm -lSDL
