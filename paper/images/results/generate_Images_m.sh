#!/bin/bash

pdflatex assembly_versus_m.tex
pdflatex error_versus_m.tex
pdflatex factor_versus_m.tex
pdflatex solve_versus_m.tex

rm *_m.aux
rm *_m.log