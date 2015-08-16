all: release

release:
	gcc findSilence.c -o findsilence -Wall -lm -lSDL

debug:
	gcc findSilence.c -o findsilence -Wall -D DEBUG \
	`pkg-config --cflags --libs gtk+-3.0` -lm -lSDL
