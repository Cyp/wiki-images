all: wiener.gif

wiener: wiener.cpp
	g++ -O2 -g -o $@ $< -Wall -Wextra -lm

wiener0099.png: wiener
	./$<

wiener.gif: wiener0099.png
	convert wiener*.png $@

clean:
	rm wiener wiener*.png
