#CXX=g++-4.3.2


PNGS = \
Orbit.png

GIFS = \
OrbitAnim.gif


all: pngs gifs
pngs: $(PNGS)
gifs: $(GIFS)

orbit: orbit.cpp ../headers/image.h
	$(CXX) -O2 -g -o $@ $< -Wall -Wextra
Orbit.ppm: orbit
	./orbit

%.jpg: %.ppm
	convert $< $@

%.png: %.ppm
	convert $< $@

OrbitAnim.gif: Orbit.ppm
	convert OrbitAnim*.ppm $@
