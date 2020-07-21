.PHONY: all clean

all: tc.png

clean:
	-rm *.pdf *.png *.log *.aux

tc.png: tc.tex cos.pdf tri.pdf rec.pdf no1.pdf no2.pdf no3.pdf
	pdflatex -shell-escape $<

%.pdf: %.tex
	pdflatex $<
