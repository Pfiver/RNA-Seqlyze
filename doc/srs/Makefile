#!/usr/bin/make

VERSION = 0.3

DIRT = $(wildcard *.out *.toc *.aux *.ps *.eps *.log *.lof *.bib *.bbl *.blg *.dvi *.tdo *~)

SRS-$(VERSION).pdf: main.tex $(wildcard *.tex *.cls)
	xelatex -halt-on-error $<
	xelatex -interaction batchmode $<
	mv main.pdf $@
	touch $@ # nudge evince
	make clean

clean:
ifneq ($(DIRT),)
	rm $(DIRT)
endif
