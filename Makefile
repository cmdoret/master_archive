include config.mk

.PHONY : all
all : $(POP)

# Running alignment with BWA
#$(MAP) : $(PROC)
#	rm -rf $@
#	mkdir -p $@
#	sed -i "s^\(main_dir=\).*^\1$(MAIN)^g" $(BWA-SRC)
#	sed -i "s/\(MM=\)[0-9]*/\1$(MM)/g" $(BWA-SRC)
#	sed -i "s/\(ALG=\)[a-z]*/\1$(ALG)/g" $(BWA-SRC)
#	sed -i "s/\(K=\)[0-9]*/\1$(K)/g" $(BWA-SRC)
#	sed -i "s/\(W=\)[0-9]*/\1$(W)/g" $(BWA-SRC)
#	bsub -K <./$(BWA-SRC)

# Running Pstacks
$(PSTACK) : $(MAP)
	bash $(P-SRC) $< $(M)
	

# Running Cstacks
$(CSTACK) : $(PSTACK)
	rm -fr $@;
	mkdir -p $@;
	sed -i "s^\(wd=\).*^\1$(MAIN)/data^g" $(C-SRC)
	sed -i "s/\(MM=\)[0-9]*/\1$(LM)/g" $(C-SRC)
	sed -i "s/^\(M=\)[0-9]*/\1$(M)/g" $(C-SRC)
	bsub -K < $(C-SRC)

# Running Sstacks
$(SSTACK) : $(CSTACK)
	sed -i "s/^\(M=\)[0-9]*/\1$(M)/g" $(S-SRC)
	bash $(S-SRC) $<

# Running populations on each family
$(POP) : $(SSTACK) $(POP-SRC)
	rm -rf $@
	mkdir -p $@
	rm -rf $(DAT)/logs/populations  # Erasing logs from previous run
	mkdir -p $(DAT)/logs/populations
	bash $(GR-SRC) $(PSTACK) $(CSTACK) $(SSTACK)
	sed -i "s^\(od=\).*^\1$(POP)^g" $(POP-SRC)
	sed -i "s/\(R=\).*/\10\.$(R)/g" $(POP-SRC)
	sed -i "s/\(D=\).*/\1$(D)/g" $(POP-SRC)
	# Changing parameters directly in the file
	bsub -K <$(POP-SRC)
	# Submitting job (-K will hang pipeline until end of job)

# Association mapping
# 1: convert vcf to ped
# 2: import data in genABEL for GWAS
#$(ASSOC) : $(POP) $(DAT)/haploid_males
#	mkdir -p $(ASSOC)
#	bash $(VCFPED) $(POP)/*.vcf
#	Rscript $(ASSOC-SRC)
	

# Rule for building lab book figures, tables and compiling Latex script
# Needs the all main steps to be run first
.PHONY : lab_book
lab_book : $(LAB) $(MISC)
	rm  -f $(LAB)/*.log $(LAB)/*.synctex* $(LAB)/*.aux $(LAB)/*.out
	# Cleaning temporary LaTeX filesjobs
	
	Rscript src/misc/assembly_stats.R $(REF)
	# Producing table of genome assembly statistics
	python2 $(MISC)/map_param.py
	# Plotting mapping statistics
	bash $(MISC)/parse_pstacks.sh
	# Parsing pstacks output into summary table
	bash $(MISC)/parse_cstacks.sh
	# Same for cstacks
	Rscript src/misc/SNP_stats.R data/SNP_lists/m2_hom_mother.txt hom_mother.txt
	Rscript src/misc/SNP_stats.R data/SNP_lists/m2_raw_CSD_candidates.txt CSD_like.txt
	# Summarizing number of SNPs
	texi2pdf -b $(LAB)/lab_book.tex -c
	mv lab_book.pdf $(LAB)
	# Compiling LaTeX report and moving it to appropriate folder


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
	bash $(MISC)/parse_VCF.sh $(POP)
	# Parsing VCF file from populations output
	python2 src/ploidy/haplo_males.py $(VCFSUM)
	# Building list of haploid males
	mkdir -p data/ploidy/plots/density
	mkdir -p data/ploidy/plots/barplots
	Rscript src/ploidy/comp_thresh.R data/ploidy/thresholds/
	# Visualize different thresholds with resulting ploidies
	mkdir -p reports/lab_book/assoc_explo_fam
	mkdir -p data/SNP_lists
	# Creating folder to store new plots if necessary
	python2 src/misc/explo_assoc.py data/ploidy/vcftools/ data/ploidy/thresholds/m2
	# Plotting exploratory results for het. at each SNP
	mkdir -p reports/lab_book/ploidy_per_fam
	# for t in data/ploidy/thresholds/*; do Rscript src/ploidy/prop_offspring.R $$t;done
	# Proportion of offspring type per family

# Saving an archive folder with all the data and parameters used.
# Not tested, probably very slow and memory consuming. Careful with this.
.PHONY : archive
archive:
	mkdir -p archive/ALG$(ALG)MM$(MM)K$(K)W$(W)_M$(M)_LM$(LM)_R$(R)D$(D)
	cp -r * archive/ALG$(ALG)MM$(MM)K$(K)W$(W)_M$(M)_LM$(LM)_R$(R)D$(D)
