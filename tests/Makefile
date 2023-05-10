all: 450k-demo.html epic-demo.html epic2-demo.html 450k-and-epic.html complex-demo.html ewas.html cnv.html random.html minfi.html cord.html GSE86831.html gds.html

%.html: %.rmd
	/usr/bin/Rscript knit.r $<

clean: 
	rm -rf 450k-demo \
		epic-demo \
		epic2-demo \
		450k-and-epic \
		complex-demo \
		GSE86831-demo \
		ewas \
		cord \
		cnv \
		random \
		minfi \
		gds \
		figure cache

## To check for any errors:
##  grep -e Error *.md */*.md

## To copy outputs to OUTPUT_DIR:
##  cp -rv random random.{md,html} cnv cnv.{md,html} 450k-and-epic 450k-and-epic.{md,html} complex-demo complex-demo.{md,html} 450k-demo 450k-demo.{md,html} epic-demo epic-demo.{md,html} epic2-demo epic2-demo.{md,html} ewas ewas.{md,html} cord cord.{md,html} minfi minfi.{md,html} figure OUTPUT_DIR
