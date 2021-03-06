include config.mk

.PHONY : all
all : $(POP)

# Running alignment with BWA
#$(MAP) : $(PROC)
#	rm -rf $@
#	mkdir -p $@
#	sed -i'' "s/\(MM=\)[0-9]*/\1$(MM)/g" $(BWA-SRC)
#	sed -i'' "s/\(ALG=\)[a-z]*/\1$(ALG)/g" $(BWA-SRC)
#	sed -i'' "s/\(K=\)[0-9]*/\1$(K)/g" $(BWA-SRC)
#	sed -i'' "s/\(W=\)[0-9]*/\1$(W)/g" $(BWA-SRC)
#	bsub -J "map_gen" -K <./$(BWA-SRC)

# Running Pstacks
#$(PSTACK) : $(MAP)
#	bash $(P-SRC) $< $(M)


# Running Cstacks
#$(CSTACK): $(PSTACK)
#	rm -fr $@;
#	mkdir -p $@;
#	sed -i'' "s^\(wd=\).*^\1$(DAT)^g" $(C-SRC)
#	sed -i'' "s/\(MM=\)[0-9]*/\1$(LM)/g" $(C-SRC)
#	sed -i'' "s/^\(M=\)[0-9]*/\1$(M)/g" $(C-SRC)
#	bsub -K < $(C-SRC)

# Running Sstacks
$(SSTACK) : $(CSTACK)
	sed -i'' "s/^\(M=\)[0-9]*/\1$(M)/g" $(S-SRC)
	sed -i'' "s^\(wd=\).*^\1$(DAT)^g" $(S-SRC)
	bash $(S-SRC) $<
	bash $(GR-SRC) $(PSTACK) $(CSTACK) $(SSTACK) $(GRFAM)

# Running populations on each family
$(POP) : $(SSTACK) $(POP-SRC)
	rm -rf $@
	mkdir -p $@
	# Erasing logs from previous run
	rm -rf $(DAT)/logs/populations
	mkdir -p $(DAT)/logs/populations
	sed -i'' "s^\(od=\).*^\1$(POP)^g" $(POP-SRC)
	sed -i'' "s/\(R=\).*/\10\.$(R)/g" $(POP-SRC)
	sed -i'' "s/\(D=\).*/\1$(D)/g" $(POP-SRC)
	sed -i'' "s^\(thresh=\).*^\1$(THRESH)^g" $(POP-SRC)
	sed -i'' "s/\(group=\).*/\1$(GRFAM)/g" $(POP-SRC)
	# Changing parameters directly in the file
	bsub -K < $(POP-SRC)
	# Submitting job (-K will hang pipeline until end of job)

# Centromere identification
$(CENTRO) :
	mkdir -p $(CENTRO)/plots
	# Processing "genomic" output from populations to get fixed and variant sites
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) $(ASSOC) $(GRFAM) --pool_output
	Rscript $(ASSOC-SRC)/chrom_types.R $(ASSOC)/grouped_outpool_prophom.tsv $(CENTRO)

# Processing genomic output for association mapping
$(ASSOC)/grouped_prophom.tsv :
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) $(ASSOC) $(GRFAM)

# Association mapping
.PHONY : assoc_mapping
assoc_mapping : $(CENTRO) $(ASSOC)/grouped_prophom.tsv
	mkdir -p $(ASSOC)/plots
	mkdir -p $(ASSOC)/hits
	Rscript $(ASSOC-SRC)/group_mothers.R $(NCSD) $(THRESH) $(ASSOC)/mother_groups.tsv
	Rscript $(ASSOC-SRC)/CSD_scan.R $(ASSOC)/grouped_prophom.tsv $(REF-ANN) $(ASSOC) 0.85
	mkdir -p $(ASSOC)/case_control
	Rscript $(ASSOC-SRC)/case_control.R $(ASSOC)/mother_groups.tsv $(ASSOC)/grouped_prophom.tsv $(ASSOC)/case_control/


# This rule is used to split haploid and diploid males.
# This has already been done under stringent parameters (d=25)
.PHONY : ploidy
ploidy:
	mkdir -p $(DAT)/ploidy/thresholds
	bash $(MISC)/parse_VCF.sh $(POP) $(GRFAM)
	# Parsing VCF file from populations output
	python2 src/ploidy/haplo_males.py $(VCFSUM) $(THRESH)
	# Building list of haploid males
	mkdir -p reports/lab_book/assoc_explo_fam
	mkdir -p data/SNP_lists
	# Creating folder to store new plots if necessary
	#python2 src/misc/explo_assoc.py data/ploidy/vcftools/ $(THRESH) $(GRFAM)
	# Plotting exploratory results for het. at each SNP
	mkdir -p reports/lab_book/ploidy_per_fam
	Rscript src/ploidy/prop_offspring.R $(THRESH)
	# Proportion of offspring type per family


# Rule for building lab book figures, tables and compiling Latex script
# Needs the all main steps to be run first
.PHONY : lab_book
lab_book : $(LAB) $(MISC)
	rm  -f $(LAB)/*.log $(LAB)/*.synctex* $(LAB)/*.aux $(LAB)/*.out
	# Cleaning temporary LaTeX filesjobs
	Rscript src/misc/assembly_stats.R $(REF)
	# Producing table of genome assembly statistics
	#bash src/mapping/parse_summaries.sh
	#python2 $(MISC)/map_param.py
	# Plotting mapping statistics
	#bash $(MISC)/parse_pstacks.sh
	# Parsing pstacks output into summary table
	#bash $(MISC)/parse_cstacks.sh
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


# Saving an archive folder with all the data and parameters used.
# Not tested, probably very slow and memory consuming. Careful with this.
.PHONY : archive
archive:
	mkdir -p archive/ALG$(ALG)MM$(MM)K$(K)W$(W)_M$(M)_LM$(LM)_R$(R)D$(D)
	cp -r * archive/ALG$(ALG)MM$(MM)K$(K)W$(W)_M$(M)_LM$(LM)_R$(R)D$(D)
