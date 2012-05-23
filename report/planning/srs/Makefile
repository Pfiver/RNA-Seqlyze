#!/usr/bin/make

VERSION = 0.1

DIRT = $(wildcard *.out *.toc *.aux *.ps *.eps *.log *.lof *.bib *.bbl *.blg *.dvi *.tdo *~)

SRS-$(VERSION).pdf: main.tex $(wildcard *.tex *.cls)
	pdflatex -halt-on-error -draftmode $< < /dev/null
	pdflatex -interaction batchmode $<
	mv main.pdf $@
	touch $@ # nudge evince
	make clean

clean:
	rm $(DIRT)
