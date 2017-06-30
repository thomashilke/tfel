PDFLATEX = pdflatex

STDOUT_REDIRECT = /dev/null
#STDOUT_REDIRECT = &1

.DEFAULT_GOAL = all
.PHONY = all pdf clean

FIGURES = figures/quadrature-convergence.tex
figures/quadrature-convergence.tex: figures/quadrature-convergence.plot\
				    data/quadrature_1d_convergence.dat

$(FIGURES): figures/%.tex:
	@echo [GNUPLOT]
	@cd figures && gnuplot $*.plot


SOURCES = src/document.tex src/content.tex

all: pdf

pdf: pdf/document.pdf

pdf/document.pdf:  build/document.pdf
	@echo "[CP]" $@
	@cp $< $@

build/document.pdf: $(SOURCES) $(FIGURES)
	@echo "[PDFLATEX]" $@
	@cd $(dir $<) && $(PDFLATEX) -shell-escape -output-directory=../$(dir $@) -halt-on-error -jobname=document $(patsubst src/%.tex,%.tex,$<) >$(STDOUT_REDIRECT)
	@cd $(dir $<) && $(PDFLATEX) -shell-escape -output-directory=../$(dir $@) -halt-on-error -jobname=document $(patsubst src/%.tex,%.tex,$<) >$(STDOUT_REDIRECT)
	@cd $(dir $<) && $(PDFLATEX) -shell-escape -output-directory=../$(dir $@) -halt-on-error -jobname=document $(patsubst src/%.tex,%.tex,$<) >$(STDOUT_REDIRECT)

clean:
	@rm -rf build/*
	@find . -name "*~" -delete


