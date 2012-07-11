#!/usr/bin/make

VERSION = 0.1

DIRT = $(wildcard *.out *.toc *.aux *.ps *.eps *.log *.lof *.bbl *.blg *.dvi *.tdo *-blx.bib *~)

SVGS = $(wildcard img/*.svg)

BachelorThesis-$(VERSION).pdf: main.tex $(wildcard *.tex *.cls *.bib) $(patsubst %.svg,%.pdf,$(SVGS))
	xelatex -halt-on-error $<
	bibtex $(patsubst %.tex,%,$<)
	xelatex -interaction batchmode $<
	xelatex -interaction batchmode $<
	mv main.pdf $@
	touch $@ # nudge evince
	make clean

sourcecode.tex:	../apidoc/_build/latex/RNA-Seqlyze.tex
	sed '/^\\subsubsection.*rnaseqlyze/,$$!d' $< | head -n-5 > $@

img/%.pdf: $(SVGS)
	convert $< $@

clean:
ifneq ($(DIRT),)
	rm $(DIRT)
endif
