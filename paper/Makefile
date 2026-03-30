TEX      := paper
FIGURES  := $(wildcard figures/*.pdf)

# PDF via pdflatex + bibtex (triple pass for references)
.PHONY: pdf html clean

pdf: $(TEX).pdf

$(TEX).pdf: $(TEX).tex $(TEX).bib $(FIGURES)
	pdflatex -interaction=nonstopmode $(TEX).tex
	bibtex $(TEX)
	pdflatex -interaction=nonstopmode $(TEX).tex
	pdflatex -interaction=nonstopmode $(TEX).tex

# HTML via tex2html (LaTeXML backend, academic theme)
html: $(TEX).html

$(TEX).html: $(TEX).tex $(TEX).bib $(FIGURES)
	tex2html -t academic -o $(TEX).html $(TEX).tex

clean:
	rm -f $(TEX).aux $(TEX).bbl $(TEX).blg $(TEX).log $(TEX).out
	rm -f $(TEX).pdf
	rm -rf $(TEX).html
	rm -f texput.log
