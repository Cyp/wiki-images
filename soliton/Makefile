CXX=g++-4.3.2


#JPGS = \
#FeynmanPropagatorWithMass0.2.jpg \
#FeynmanPropagatorWithMass2.jpg \
#FeynmanPropagatorWithMass20.jpg \
#FeynmanPropagatorWithMass200.jpg \
#PropagatorColours.jpg

PNGS = \
DoubleWellSoliton.png

GIFS = \
DoubleWellSolitonAntisoliton.gif

all: pngs gifs #jpgs
$(JPGS): FeynmanPropagatorWithMass2.ppm
#jpgs: $(JPGS)
pngs: $(PNGS)
gifs: $(GIFS)

#-fopenmp
#field: field.cpp
#	$(CXX) -O2 -s -o $@ $< -Wall -Wextra

#FeynmanPropagatorWithMass2.ppm: field
#	./field

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
