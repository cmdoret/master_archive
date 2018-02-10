include config.mk

ifdef $(LOCAL)
	LOCAL='$(LOCAL)';
endif

# reference-based STACKS pipeline for LSF cluster
ref_lsf:
	make -f src/pipelines/Makeref.lsf

# reference based without LSF
.PHONY : ref_nix
ref_nix:
	make -f src/pipelines/Makeref.nix


# Centromere identification
$(CENTRO) :
	mkdir -p $(CENTRO)/plots
	# Processing "genomic" output from populations to get fixed and variant sites
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) $(ASSOC) $(GRFAM) \
	                                        --pool_output
	# Inferring centromere position based on recombination raes along chromosomes
	Rscript $(ASSOC-SRC)/chrom_types.R $(ASSOC)/grouped_outpool_prophom.tsv \
	                                   $(CENTRO)

# Processing genomic output for association mapping
$(ASSOC)/grouped_prophom.tsv :
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) $(ASSOC) $(GRFAM) \
	                                        --keep_all --pool_output

# Association mapping
.PHONY : assoc_mapping
assoc_mapping : $(CENTRO) $(ASSOC)/grouped_prophom.tsv
  # Creating folder to store association mapping results
	mkdir -p $(ASSOC)/plots
	mkdir -p $(ASSOC)/hits
	# Genome-wide association mapping
	mkdir -p $(ASSOC)/case_control
	Rscript $(ASSOC-SRC)/case_control.R \
					$(ASSOC)/grouped_outpool_keep_prophom.tsv \
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
	                                  $(THRESH) \
									                  --ploidy_thresh 0.90
  # Processing populations genomic output (excluding hom/missing SNPs in
	# mothers from their families)
	python2 $(ASSOC-SRC)/process_genomic.py $(POP) \
	                                        $(DAT)/ploidy/ \
																					$(GRFAM)
  # Blacklisting loci that are heterozygous in haploid males
	python2 src/ploidy/blacklist_haploloci.py $(DAT)/ploidy/ \
	                               $(BLACK) \
																 $(THRESH) \
																 $(GRFAM)

	mkdir -p data/SNP_lists

# Assembling transcripts and measuring coverage along genome
$(RNA)/assembled/ :
	bash $(RNA-SRC) -a $(BAM) -r $(OLD-REF) -o $@

# Preparing input file for collinearity analysis. Run MCScanX on files manually
.PHONY : mult_align
mult_align : $(RNA)/assembled/
	#LOCAL='$(LOCAL)'
	rm -rf $(MCSX-IN)/
	mkdir -p $(MCSX-IN)/
	bash $(MCSX-SRC) -g $(RNA)/assembled/transcripts.gtf \
	                 -o $(MCSX-IN) \
									 -r $(REF) \
									 -c $(CORRESP) \
									 $${LOCAL:+-l}
	bash src/circos_conf/circos_input_gen.sh $(REF)

# Rule for building lab book
.PHONY : lab_book
lab_book : $(LAB) $(MISC)
	rm  -f $(LAB)/*.log $(LAB)/*.synctex* $(LAB)/*.aux $(LAB)/*.out
	# Cleaning temporary LaTeX files
	Rscript src/misc/assembly_stats.R $(REF)
	# Producing table of genome assembly statistics
	texi2pdf -b $(LAB)/lab_book.tex -c
	mv lab_book.pdf $(LAB)
	# Compiling LaTeX report and moving it to appropriate folder

.PHONY : wgs_wild_mothers
wgs_wild_mothers :
	#bash src/wgs_wild_mothers/bwa_wgs.sh --workdir $(WGS) --ref $(REF) --out
	bash src/wgs_wild_mothers/snps_wgs.sh --workdir $(WGS) --ref $(REF) --winsize 100

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
