.PHONY:	all postscript patch_index ascii figs pdffigs pdf tpdf clean veryclean world

DOC = manual

all:	$(DOC).tex
	make ascii postscript

world:	$(DOC).tex install-doc
	make veryclean
	make tpdf
	make clean
	make postscript
	make clean
	make ascii

install_mod_date=`ls --time-style=+'%Y %b %d' -l ../INSTALL.txt | cut -d' ' -f6,7,8`

install-doc:	../INSTALL.html ../INSTALL.xml ../INSTALL.pdf
../INSTALL.html:	../INSTALL.txt
	asciidoc -a "revdate=${install_mod_date}" ../INSTALL.txt
../INSTALL.xml:	../INSTALL.txt
	asciidoc -a "revdate=${install_mod_date}" -b docbook ../INSTALL.txt
../INSTALL.pdf:	../INSTALL.txt
	dblatex -T simple ../INSTALL.xml

TOC=\\begin{theindex}\\label{chap:index}\\ifpdf\\pdfbookmark[-1]{Index}{chap:index}\\fi

patch_index:
	sed -e 's/\\begin{theindex}/$(TOC)/' $(DOC).ind > tmp.ind
	mv -f tmp.ind $(DOC).ind

postscript:	$(DOC).tex
	pdftops $(DOC).pdf

old-postscript:	$(DOC).tex
	latex -output-format=dvi $(DOC).tex
	latex -output-format=dvi $(DOC).tex
	makeindex $(DOC).idx
	make patch_index
	latex -output-format=dvi $(DOC).tex
	latex -output-format=dvi $(DOC).tex
	dvips -o $(DOC).ps $(DOC)

ascii:	$(DOC).tex
	jed -batch -l ./extract_help.sl

scripts = contin.sl crude_dem.sl custom1.sl data.sl elev.sl \
          ion_curves.sl line_id.sl profile_fit.sl xspec.sl \
	  residuals.sl contour.sl

figs:
	-mkdir figures
	cd figures; for s in $(scripts) ; do \
	   echo "Doing $$s"; isis --batch ../scripts/$$s; \
	done

fignames = contin data_flx line_id crude_dem elev profile_fit \
           custom1 ion_bal residuals data_cts ion_frac xspec contour

pdffigs:
	cd figures; for s in $(fignames) ; do \
	   ps2epsi $$s.ps ; epstopdf $$s.epsi ;\
	done

pdf:
	pdflatex $(DOC)
	pdflatex $(DOC)
	makeindex $(DOC).idx
	make patch_index
	pdflatex $(DOC)
	pdflatex $(DOC)

tpdf:	$(DOC).pdf
	thumbpdf $(DOC)
	pdflatex $(DOC)

$(DOC).pdf:
	make pdf

clean:
	/bin/rm -f $(DOC).aux $(DOC).dvi $(DOC).idx $(DOC).ilg $(DOC).out \
                $(DOC).ind $(DOC).lof $(DOC).log $(DOC).lot $(DOC).toc \
	        $(DOC).tpt thumb*.??? index.log

veryclean:
	make clean
	/bin/rm -f $(DOC).ps $(DOC).pdf help.txt
	/bin/rm -f ../INSTALL.html ../INSTALL.pdf ../INSTALL.xml

