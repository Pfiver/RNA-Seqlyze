#!/usr/bin/make

VERSION = 0.1

DIRT = $(wildcard *.out *.toc *.aux *.ps *.eps *.log *.lof *.bib *.bbl *.blg *.dvi *.tdo *~)

BachelorThesis-$(VERSION).pdf: main.tex $(wildcard *.tex *.cls)
	xelatex -halt-on-error $<
	xelatex -interaction batchmode $<
	mv main.pdf $@
	touch $@ # nudge evince
	make clean

sourcecode.tex:	../apidoc/_build/latex/RNA-Seqlyze.tex
	sed '/^\\subsubsection.*rnaseqlyze/,$$!d' $< | head -n-5 > $@

clean:
ifneq ($(DIRT),)
	rm $(DIRT)
endif
