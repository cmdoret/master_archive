include config.mk

.PHONY : all
all : $(POP)

# Running alignment with BWA
#$(MAP) : $(PROC)
# 	rm -rf $@
# 	mkdir -p $@
# 	sed -i "s^\(main_dir=\).*^\1$(MAIN)^g" $(BWA-SRC)
# 	sed -i "s/\(MM=\)[0-9]*/\1$(MM)/g" $(BWA-SRC)
# 	sed -i "s/\(ALG=\)[a-z]*/\1$(ALG)/g" $(BWA-SRC)
# 	sed -i "s/\(K=\)[0-9]*/\1$(K)/g" $(BWA-SRC)
# 	sed -i "s/\(W=\)[0-9]*/\1$(W)/g" $(BWA-SRC)
# 	bsub -K <./$(BWA-SRC)

# Running Pstacks
# $(PSTACK) : $(MAP)
# 	bash $(P-SRC) $< $(M)
	

# Running Cstacks
# $(CSTACK) : $(PSTACK)
# 	rm -fr $@;
# 	mkdir -p $@;
# 	sed -i "s^\(wd=\).*^\1$(MAIN)/data^g" $(C-SRC)
# 	sed -i "s/\(MM=\)[0-9]*/\1$(LM)/g" $(C-SRC)
# 	sed -i "s/^\(M=\)[0-9]*/\1$(M)/g" $(C-SRC)
# 	bsub -K < $(C-SRC)

# Running Sstacks
$(SSTACK) : $(CSTACK)
	sed -i "s/^\(M=\)[0-9]*/\1$(M)/g" $(S-SRC)
	bash $(S-SRC) $<

# Running populations
$(POP) : $(SSTACK) $(POP-SRC)
	rm -rf $@
	mkdir -p $@
	rm -rf $(DAT)/logs/populations
	mkdir -p $(DAT)/logs/populations
	bash $(GR-SRC) $(PSTACK) $(CSTACK) $(SSTACK)
	sed -i "s^\(od=\).*^\1$(POP)^g" $(POP-SRC)
	sed -i "s/\(R=\).*/\10\.$(R)/g" $(POP-SRC)
	sed -i "s/\(D=\).*/\1$(D)/g" $(POP-SRC)
	bsub -K <$(POP-SRC)
	mv $(SSTACK)/batch* $(POP)

# Association mapping
# 1: convert vcf to ped
# 2: import data in genABEL for GWAS
$(ASSOC) : $(POP) $(DAT)/haploid_males
	mkdir -p $(ASSOC)
	bash $(VCFPED) $(POP)/*.vcf
	Rscript $(ASSOC-SRC)
	

# Rule for building lab book figures, tables and compiling Latex script
# Needs the all main steps to be run first
.PHONY : lab_book
lab_book : $(LAB) $(MISC)
	rm  -f $(LAB)/*.log $(LAB)/*.synctex* $(LAB)/*.aux $(LAB)/*.out
	Rscript src/misc/assembly_stats.R $(REF)
	python2 $(MISC)/map_param.py
	bash $(MISC)/parse_pstacks.sh
	bash $(MISC)/parse_cstacks.sh
	bash $(MISC)/parse_VCF.sh
	Rscript $(MISC)/plot_VCF.R $(MISC)
	mkdir -p reports/lab_book/assoc_explo
	for t in data/ploidy/thresholds/*; do python2 src/misc/explo_assoc.py $(POP)/*haplotypes.tsv $$t;done
	mkdir -p reports/lab_book/ploidy_per_fam
	for t in data/ploidy/thresholds/*; do Rscript src/ploidy/prop_offspring.R $$t;done
	texi2pdf -b $(LAB)/lab_book.tex -c
	mv lab_book.pdf $(LAB)


# When it gets too messy
.PHONY : clean
clean :
	rm -f *STDERR*
	rm -f *STDOUT*
	rm -f demulti*
	rm -rf bam
	rm -rf bsub_scripts


# This rule is used to split haploid and diploid males.
# This has already been done under stringent parameters (d=25)
.PHONY : ploidy
ploidy:
	mkdir -p $(DAT)/ploidy/thresholds
	bash $(MISC)/parse_VCF.sh
	python2 src/ploidy/haplo_males.py $(VCFSUM)
	mkdir -p data/ploidy/plots
	Rscript src/ploidy/comp_thresh.R data/ploidy/thresholds/
