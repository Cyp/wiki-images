CXX=g++-4.3.2


JPGS = \
FeynmanPropagatorWithMass0.2.jpg \
FeynmanPropagatorWithMass2.jpg \
FeynmanPropagatorWithMass20.jpg \
FeynmanPropagatorWithMass200.jpg \
PropagatorColours.jpg


all: jpgs
$(JPGS): FeynmanPropagatorWithMass2.ppm
jpgs: $(JPGS)

#-fopenmp
field: field.cpp
	$(CXX) -O2 -s -o $@ $< -Wall -Wextra

FeynmanPropagatorWithMass2.ppm: field
	./field

%.jpg: %.ppm
	convert $< $@
