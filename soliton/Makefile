CXX=g++-4.3.2


PNGS = \
DoubleWellSoliton.png

GIFS = \
DoubleWellSolitonAntisoliton.gif

all: pngs gifs
$(JPGS): FeynmanPropagatorWithMass2.ppm
pngs: $(PNGS)
gifs: $(GIFS)

soliton: soliton.cpp
	$(CXX) -O2 -s -o $@ $< -Wall -Wextra
soliton.ppm: soliton
	./soliton

%.jpg: %.ppm
	convert $< $@

DoubleWellSoliton.png: soliton.ppm
	convert $< $@

DoubleWellSolitonAntisoliton.gif: soliton.ppm
	convert disoliton*.ppm $@
