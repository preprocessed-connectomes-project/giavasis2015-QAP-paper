
all: qap.tex qap.bib Makefile
	pdflatex qap.tex
	bibtex qap
	pdflatex qap.tex
	pdflatex qap.tex

clean:
	rm -f *.aux *.bbl *.log *.dvi
