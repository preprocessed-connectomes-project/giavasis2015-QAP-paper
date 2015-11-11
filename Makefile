
all: main.tex qap.bib Makefile
	pdflatex main.tex
	bibtex main
	pdflatex main.tex
	pdflatex main.tex

clean:
	rm -f *.aux *.bbl
