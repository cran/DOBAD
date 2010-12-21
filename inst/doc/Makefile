# HORRIBLY DANGEROUS: this removes ALL TEX FILES! So they don't ship with the package.
# [I have now removed  the delete-the-latex-files line from 'all'; the  tex files might as well ship with the package]
#rm -rf ${TEXSOURCE}
#rm -rf .tex *.ps *.dvi
#	pdflatex ${TEXSOURCE}	##THIS DOENST WORK FOR SOME REASON
# ${SOURCE:.tex=.pdf}

.SUFFIXES: .Rnw .tex .dvi .ps .pdf
PDFS = BD_EM.pdf BD_EM_wImmig.pdf BD_MCMC.pdf sim_condBD.pdf

all: ${PDFS}
	rm -rf *.ps *.dvi *.aux *.bbl *.blg *.log

.Rnw.tex:
	R CMD Sweave $<

.tex.dvi:
	latex $<
	bibtex $*
	latex $<
	latex $<

.dvi.ps:
	dvips -t letter $< -o

.ps.pdf: 
	ps2pdf $<

clean:
	rm -rf *.ps *.dvi *.aux *.bbl *.blg *.log
