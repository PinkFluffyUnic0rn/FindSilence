all: release

release:
	gcc findSilence.c -o findsilence -lm -lSDL

debug:
	gcc findSilence.c -o findsilence \
	`pkg-config --cflags --libs gtk+-3.0` -lm -lSDL
