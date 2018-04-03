all: 450k-demo.html epic-demo.html 450k-and-epic.html ewas.html cnv.html random.html minfi.html cord.html GSE86831.html gds.html

%.html: %.rmd
	Rscript knit.r $<

clean: 
	rm -rf 450k-demo 450k-demo.{md,html} \
		epic-demo epic-demo.{md,html} \
		450k-and-epic 450k-and-epic.{md,html}  \
		GSE86831-demo GSE86831.{md,html} \
		ewas ewas.{md,html} \
		cord cord.{md,html} \
		cnv cnv.{md,html} \
		random random.{md,html} \
		minfi minfi.{md,html} \
		gds gds.{md,html} \
		figure cache

## To check for any errors:
##  grep -e Error *.md */*.md

## To copy outputs to OUTPUT_DIR:
##  cp -rv random random.{md,html} cnv cnv.{md,html} 450k-and-epic 450k-and-epic.{md,html} 450k-demo 450k-demo.{md,html} epic-demo epic-demo.{md,html} ewas ewas.{md,html} cord cord.{md,html} minfi minfi.{md,html} figure OUTPUT_DIR
