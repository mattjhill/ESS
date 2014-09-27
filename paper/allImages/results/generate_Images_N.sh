#!/bin/bash

pdflatex assembly_versus_N.tex
pdflatex error_versus_N.tex
pdflatex factor_versus_N.tex
pdflatex solve_versus_N.tex

rm *_N.aux
rm *_N.log