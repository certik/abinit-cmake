#
# Makefile for the theory documents of ABINIT
#

# Names of the documents
documents = \
  1WF.pdf \
  geometry.pdf \
  gwa.pdf \
  noncol.pdf \
  set_WFs.pdf \
  vanderbilt-anaddb-notes.pdf

# --------------------------------------------------------------------------- #

#
# Rules to build the documents
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

all_targets all: $(documents)

clean:
	rm -f *.tmp *.blg *.ilg *.log *.aux *.nav *.out *.snm *.toc *.vrb

mostlyclean: clean
	rm -f *.dvi

distclean: mostlyclean
	rm -f $(documents)
