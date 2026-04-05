TEX      := paper
SUPP     := supplementary
FIGURES  := $(wildcard figures/*.pdf)

.PHONY: pdf supp all html clean

all: pdf supp

# Main paper via pdflatex + bibtex (triple pass for references)
pdf: $(TEX).pdf

$(TEX).pdf: $(TEX).tex $(TEX).bib $(FIGURES)
	pdflatex -interaction=nonstopmode $(TEX).tex
	bibtex $(TEX)
	pdflatex -interaction=nonstopmode $(TEX).tex
	pdflatex -interaction=nonstopmode $(TEX).tex

# Supplementary material
supp: $(SUPP).pdf

$(SUPP).pdf: $(SUPP).tex
	pdflatex -interaction=nonstopmode $(SUPP).tex
	pdflatex -interaction=nonstopmode $(SUPP).tex

# HTML via tex2html (LaTeXML backend, academic theme)
html: $(TEX).html

$(TEX).html: $(TEX).tex $(TEX).bib $(FIGURES)
	tex2html -t academic -o $(TEX).html $(TEX).tex

clean:
	rm -f $(TEX).aux $(TEX).bbl $(TEX).blg $(TEX).log $(TEX).out
	rm -f $(TEX).pdf
	rm -f $(SUPP).aux $(SUPP).log $(SUPP).out $(SUPP).pdf
	rm -rf $(TEX).html
	rm -f texput.log
