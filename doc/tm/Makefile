# -*- sh -*-
#
# To create the SGML files, you will need to install the tm-utils
# package.  See http://www.jedsoft.org/ for more information.
#
TMEXPAND = /aluche/d1/web/tm-dist/bin/tmexpand
SL2TM = sl2tm
MACRODIR = /aluche/d1/web/tm-dist/macros

MODULE = cfitsio
HLPFUNS_TM = rtl/cfitsiofuns.tm fitsfuns.tm fitswcsfuns.tm
AUTOGEN_TM = fitswcsfuns.tm fitsfuns.tm
MODULE_DEPS = rtl/cfitsiofuns.tm $(AUTOGEN_TM)

TXT_FILES = $(MODULE).txt
SGML_FILES = $(MODULE).sgml
HTML_FILES = $(MODULE).html
TEX_FILES = $(MODULE).tex
PS_FILES = $(MODULE).ps 
PDF_FILES = $(MODULE).pdf
HLP_FILE = $(MODULE).hlp

HLP_TM_FILE = helpfile.tm

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
HELPDIR = ../help

SUBDIRS = $(TEXTDIR) $(HTMLDIR) $(PSDIR) $(SGMLDIR) $(PDFDIR) $(HELPDIR)
SRCDIR = `pwd`

all: $(SGML_FILES) $(HTML_FILES) $(TEX_FILES) $(TXT_FILES)
fitsfuns.tm: ../../src/fits.sl
	$(SL2TM) ../../src/fits.sl > fitsfuns.tm
fitswcsfuns.tm: ../../src/fitswcs.sl
	$(SL2TM) ../../src/fitswcs.sl > fitswcsfuns.tm
text-files: $(TXT_FILES)
#----- SGML Files -----------------------------------------------------------
$(MODULE).sgml : $(MODULE).tm $(MODULE_DEPS)
	$(TMEXPAND) -I$(MACRODIR) $(MODULE).tm $(MODULE).sgml
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
help-files: $(HLP_FILE)
$(HLP_FILE): $(HLP_TM_FILE) $(HLPFUNS_TM)
	$(TMEXPAND) -I$(MACRODIR) $(HLP_TM_FILE) $(HLP_FILE)
texclean:
	-rm -f *.dvi *.log *.aux *.toc *.out
clean: texclean
	-rm -f *~ rtl/*.BAK rtl/*~ *.tmp *-error
distclean: clean
	-rm -f *.html *.ps $(HLP_FILE) $(TXT_FILES) $(TEX_FILES) $(SGML_FILES) $(PDF_FILES) $(AUTOGEN_TM)
install-txt: $(TXT_FILES)
	-mv $(TXT_FILES) ../text
install-help: $(HLP_FILE)
	-mkdir -p $(HELPDIR)
	-mv $(HLP_FILE) $(HELPDIR)
install-all: all install-help install-txt $(PS_FILES) $(PDF_FILES)
	-mkdir -p $(HTMLDIR) $(PSDIR) $(SGMLDIR) $(PDFDIR)
	-mv *.html $(HTMLDIR)
	-mv $(PS_FILES) ../ps
	-mv $(SGML_FILES) ../sgml
	-mv $(PDF_FILES) ../pdf
install: install-txt install-help
