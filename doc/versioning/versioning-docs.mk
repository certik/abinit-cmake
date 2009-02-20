#
# Makefile for the source code lecture of the ABINIT 2005 summer school
#

# Name of the documents
documents = \
 bzr-quickref.pdf \
 bzr-intro-p1.html \
 bzr-intro-p2.html \
 bzr-intro-p3.html \
 vcs-comparison.html

# --------------------------------------------------------------------------- #

#
# Rules to build the document
#

LATEX    = latex
TEXFLAGS = 
DVI2PDF  = dvipdf
D2PFLAGS = 
MARKDOWN = markdown
MDFLAGS  = 

.SUFFIXES:
.SUFFIXES: .pdf .dvi .tex .html .txt

.tex.dvi:
	$(LATEX) $(TEXFLAGS) $<
	$(LATEX) $(TEXFLAGS) $<

.dvi.pdf:
	$(DVI2PDF) $(D2PFLAGS) $<

.txt.html:
	$(MARKDOWN) $(MDFLAGS) $< > $@

# --------------------------------------------------------------------------- #

#
# Targets
#

all_targets all: $(documents)

clean:
	rm -f *.tmp *.blg *.ilg *.log *.aux *.nav *.out *.snm *.toc *.vrb

mostlyclean: clean
	rm -f bzr-quickref.dvi

distclean: mostlyclean
	rm -f $(documents)

