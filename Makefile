# latex Makefile
LATEX=/usr/texbin/latex
BIBTEX=/usr/texbin/bibtex
DVIPS=/usr/texbin/dvips
TARGETS= ms.ps ms.pdf ms-aas

all: $(TARGETS)

clean: 
	@rm -f *.aux *.bbl *.blg *.dvi *.log

ms.dvi: ms.tex
	${LATEX} ms.tex

ms.ps: ms.dvi 
	${BIBTEX} ms
	${LATEX} ms.tex
	${BIBTEX} ms
	${LATEX} ms.tex
	${DVIPS} ms.dvi -o ms.ps

ms.pdf: ms.ps
	ps2pdf ms.ps

ms-aas: 
	./nat2jour.pl -references ms 
	cp ms.bbl ms-aas.bbl
	${LATEX} ms-aas.tex
	${LATEX} ms-aas.tex
	${DVIPS} ms-aas.dvi -o ms-aas.ps
	ps2pdf ms-aas.ps
ms-arxiv:
	${LATEX} ms.tex
	${LATEX} ms.tex
	${DVIPS} ms.dvi -o ms.ps
	ps2pdf ms.ps










