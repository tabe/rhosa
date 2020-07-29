.PHONY: all clean mostlyclean

all: tc.png

clean: mostlyclean
	-rm tc.png

mostlyclean:
	-rm tc.pdf cos.pdf tri.pdf rec.pdf no1.pdf no2.pdf no3.pdf *.log *.aux

tc.png: tc.tex cos.pdf tri.pdf rec.pdf no1.pdf no2.pdf no3.pdf
	pdflatex -shell-escape $<

%.pdf: %.tex
	pdflatex $<
