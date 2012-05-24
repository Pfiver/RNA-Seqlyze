#!/usr/bin/make

VERSION = 0.1

DIRT = $(wildcard *.out *.toc *.aux *.ps *.eps *.log *.lof *.bib *.bbl *.blg *.dvi *.tdo *~)

SRS-$(VERSION).pdf: main.tex $(wildcard *.tex *.cls)
	xelatex -halt-on-error -interaction batchmode $<
	xelatex -halt-on-error $<
	mv main.pdf $@
	touch $@ # nudge evince
	make clean

clean:
	rm $(DIRT)
