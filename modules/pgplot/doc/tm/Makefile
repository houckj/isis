# -*- sh -*-
#
# To create the SGML files, you will need to install the tm-utils
# package.  See http://www.jedsoft.org/ for more information.
#
TM2SGML = /aluche/d1/web/tm-dist/bin/tmexpand
MACRODIR = /aluche/d1/web/tm-dist/macros

MODULE = maplib

TXT_FILES = $(MODULE).txt
SGML_FILES = $(MODULE).sgml
HTML_FILES = $(MODULE).html
TEX_FILES = $(MODULE).tex
PS_FILES = $(MODULE).ps 
PDF_FILES = $(MODULE).pdf

SGML2LATEX = sgml2latex -p letter -o tex
SGML2HTML = sgml2html
SGML2TXT = sgml2txt -f
LATEX = latex
PDFLATEX = pdflatex

TEXTDIR = ../text
PSDIR = ../ps
HTMLDIR = ../html
SGMLDIR = ../sgml
PDFDIR = ../pdf

SUBDIRS = $(TEXTDIR) $(HTMLDIR) $(PSDIR) $(SGMLDIR) $(PDFDIR)
SRCDIR = `pwd`

all: $(SGML_FILES) $(HTML_FILES) $(TEX_FILES) $(TXT_FILES)
text-files: $(TXT_FILES)
#----- SGML Files -----------------------------------------------------------
$(MODULE).sgml : $(MODULE).tm
	$(TM2SGML) -I$(MACRODIR) $(MODULE).tm $(MODULE).sgml
#----- HTML Files -----------------------------------------------------------
$(MODULE).html : $(MODULE).sgml
	$(SGML2HTML) $(MODULE).sgml
#----- TeX Files ------------------------------------------------------------
$(MODULE).tex : $(MODULE).sgml
	$(SGML2LATEX) $(MODULE).sgml
	jed -script ./fixtex.sl $(MODULE).tex
#----- PDF Files -----------------------------------------------------------
$(MODULE).pdf : $(MODULE).tex
	$(MAKE) texclean
	$(PDFLATEX) $(MODULE).tex
	$(PDFLATEX) $(MODULE).tex
	$(PDFLATEX) $(MODULE).tex
#----- PS Files -----------------------------------------------------------
$(MODULE).ps : $(MODULE).tex texclean
	$(LATEX) $(MODULE).tex
	$(LATEX) $(MODULE).tex
	$(LATEX) $(MODULE).tex
	dvips -o $(MODULE).ps $(MODULE).dvi
#----- Text Files -----------------------------------------------------------
$(MODULE).txt: $(MODULE).sgml
	$(SGML2TXT) $(MODULE).sgml
#----------------------------------------------------------------------------
texclean:
	-rm -f *.dvi *.log *.aux *.toc *.out
clean: texclean
	-rm -f *~ rtl/*.BAK rtl/*~ *.tmp
distclean: clean
	-rm -f *.html *.ps $(TXT_FILES) $(TEX_FILES) $(SGML_FILES) $(PDF_FILES)
install-txt-files: $(TXT_FILES)
	-mv $(TXT_FILES) ../text

install: all $(PS_FILES) $(PDF_FILES)
	-mkdir $(SUBDIRS)
	-mv $(TXT_FILES) ../text
	-mv *.html ../html
	-mv $(PS_FILES) ../ps
	-mv $(SGML_FILES) ../sgml
	-mv $(PDF_FILES) ../pdf
