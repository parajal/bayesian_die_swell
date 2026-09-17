
.SUFFIXES:

.SUFFIXES: .pdf .tex

.tex.pdf:
	pdflatex $<
	pdflatex $<

all: ../userguide.pdf

../userguide.pdf: userguide.pdf
	cp -f userguide.pdf ..

clean: 
	rm -f *.dvi *.log *.blg *.pdf *.out *.toc *.aux
