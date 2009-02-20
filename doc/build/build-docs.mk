#
# Makefile for the build system guide of ABINIT 5
#

# Name of the document
document = build-system-guide

# --------------------------------------------------------------------------- #

#
# Rules to build the document
#

LATEX    = latex
TEXFLAGS = 
DVI2PDF  = dvipdf
D2PFLAGS = 

.SUFFIXES:
.SUFFIXES: .pdf .dvi .tex

.tex.dvi:
	$(LATEX) $(TEXFLAGS) $<
	$(LATEX) $(TEXFLAGS) $<

.dvi.pdf:
	$(DVI2PDF) $(D2PFLAGS) $<

# --------------------------------------------------------------------------- #

#
# Targets
#

all_targets all: $(document).pdf

clean:
	rm -f *.tmp *.blg *.ilg *.log *.aux *.nav *.out *.snm *.toc *.vrb

mostlyclean: clean
	rm -f $(document).dvi

distclean: mostlyclean
	rm -f $(document).pdf
