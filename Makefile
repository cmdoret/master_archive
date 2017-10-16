include config.mk

# reference-based STACKS pipeline for LSF cluster
ref_lsf:
	make -f src/pipelines/Makeref.lsf

# reference based without LSF
.PHONY : ref_nix
ref_nix:
	make -f src/pipelines/Makeref.nix

# denovo with LSF
.PHONY : denovo_lsf
denovo_lsf:
	make -f src/pipelines/Makenovo.lsf

# denovo without LSF
#denovo_nix:
#	make -f src/pipelines/Makenovo.nix

# Centromere identification
$(CENTRO) :
	mkdir -p $(CENTRO)/plots
	# Processing "genomic" output from populations to get fixed and variant sites
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) $(ASSOC) $(GRFAM) --pool_output
	# Inferring centromere position based on recombination raes along chromosomes
	Rscript $(ASSOC-SRC)/chrom_types.R $(ASSOC)/grouped_outpool_prophom.tsv $(CENTRO)

# Processing genomic output for association mapping
$(ASSOC)/grouped_prophom.tsv :
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) $(ASSOC) $(GRFAM)

# Association mapping
.PHONY : assoc_mapping
assoc_mapping : $(CENTRO) $(ASSOC)/grouped_prophom.tsv
  # Creating folder to store association mapping results
	mkdir -p $(ASSOC)/plots
	mkdir -p $(ASSOC)/hits
	# Clustering families based on proportion of males among diploids
	Rscript $(ASSOC-SRC)/group_mothers.R \
	        $(NCSD) \
					$(THRESH) \
					$(ASSOC)/mother_groups.tsv
	# Scanning for CSD based on heterozygosity along genome
	Rscript $(ASSOC-SRC)/CSD_scan.R \
	        $(ASSOC)/grouped_prophom.tsv \
					$(REF-ANN) \
					$(ASSOC) 0.85
	# Genome-wide association mapping
	mkdir -p $(ASSOC)/case_control
	Rscript $(ASSOC-SRC)/case_control.R \
	        $(ASSOC)/mother_groups.tsv \
					$(ASSOC)/grouped_prophom.tsv \
					$(ASSOC)/case_control/


# This rule is used to split haploid and diploid males.
# This has already been done under stringent parameters (d=25)
.PHONY : ploidy
ploidy:
	# Parsing VCF file from populations output
	mkdir -p $(DAT)/ploidy/thresholds
	bash $(MISC)/parse_VCF.sh $(POP) \
	                          $(GRFAM)
	# Building list of haploid males
	python2 src/ploidy/haplo_males.py $(VCFSUM) \
	                                  $(THRESH)
  # Processing populations genomic output
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) \
	                                        $(DAT)/ploidy/ \
																					$(GRFAM)
  # Blacklisting loci that are heterozygous in haploid males
	python2 src/ploidy/blacklist_haploloci.py $(DAT)/ploidy/ \
	                               $(BLACK) \
																 $(THRESH) \
																 $(GRFAM)
	# Creating folder to store new plots if necessary
	mkdir -p reports/lab_book/assoc_explo_fam
	mkdir -p data/SNP_lists
	# Plotting exploratory results for het. at each SNP
	#python2 src/misc/explo_assoc.py data/ploidy/vcftools/ $(THRESH) $(GRFAM)
	mkdir -p reports/lab_book/ploidy_per_fam
	# Proportion of offspring type per family
	Rscript src/ploidy/prop_offspring.R $(THRESH)


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
