
all: qap_fn.tex qap.bib Makefile
	pdflatex qap_fn.tex
	bibtex qap_fn
	pdflatex qap_fn.tex
	pdflatex qap_fn.tex

clean:
	rm -f *.aux *.bbl *.log *.dvi *.blg *.out
