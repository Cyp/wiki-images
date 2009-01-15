CFLAGS = -Os -pipe -Wall -Wextra
LIBS = -lm

GIFS = \
Lorentz_transform_of_world_line.gif \
Galilean_transform_of_world_line.gif

all: gifs
gifs: $(GIFS)

galilean: galilean.c
	$(CC) -o $@ $< $(LIBS)

lorentz: lorentz.c
	$(CC) -o $@ $< $(LIBS)

lor0099.pgm: lorentz
	./lorentz

Lorentz_transform_of_world_line.gif: lor0099.pgm
	convert lor*.pgm $@

gal0099.pgm: galilean
	./galilean

Galilean_transform_of_world_line.gif: gal0099.pgm
	convert gal*.pgm $@
