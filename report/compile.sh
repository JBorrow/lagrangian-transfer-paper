rm -rf compiled

mkdir compiled

pdflatex -output-directory="compiled" report.tex
bibtex compiled/report.aux
pdflatex -output-directory="compiled" report.tex
pdflatex -output-directory="compiled" report.tex


