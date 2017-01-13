TEX=pdflatex
BIB=bibtex

paper: proj2.tex bibliography.bib asmems4.bst asme2e.cls
	$(TEX) proj2.tex
	$(BIB) proj2
	$(TEX) proj2.tex
	$(TEX) proj2.tex

clean:
	rm -rf *.log *.aux *.pdf *~ \#* *.bbl *.blg
