include config.mk

ifdef $(LOCAL)
	LOCAL='$(LOCAL)';
endif

##########################################
#### 1. RADseq processing with STACKS ####
##########################################

# reference-based STACKS pipeline for LSF cluster
ref_lsf:
	make -f src/pipelines/Makeref.lsf

# reference based STACKS without LSF
.PHONY : ref_nix
ref_local:
	make -f src/pipelines/Makeref.local

##############################
#### 2. PLOIDY SEPARATION ####
##############################

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

################################
#### 3. ASSOCIATION MAPPING ####
################################

# Association mapping
.PHONY : assoc_mapping
assoc_mapping : $(CENTRO) $(SIZES) $(ASSOC)/grouped_prophom.tsv
  # Creating folder to store association mapping results
	mkdir -p $(ASSOC)/plots
	mkdir -p $(ASSOC)/hits
	# Genome-wide association mapping
	mkdir -p $(ASSOC)/case_control
	Rscript $(ASSOC-SRC)/case_control.R \
					$(ASSOC)/grouped_outpool_keep_prophom.tsv \
					$(ASSOC)/case_control/


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


# Store sizes of all contigs in text file
$(SIZES):
	awk '$$0 ~ ">" {print c; c=0;printf substr($$0,2,100) "\t"; } \
			 $$0 !~ ">" {c+=length($$0);} END { print c; }' $(REF) > $@

##################################
#### 4. COLLINEARITY ANALYSIS ####
##################################

# Preparing input file for collinearity analysis. Run MCScanX on files manually
.PHONY : collinearity
mult_align : $(RNA)/assembled/ $(CORRESP)
	#LOCAL='$(LOCAL)'
	rm -rf $(MCSX-IN)/
	mkdir -p $(MCSX-IN)/
	bash $(MCSX-SRC) -g $(RNA)/assembled/transcripts.gtf \
	                 -o $(MCSX-IN) \
									 -r $(REF) \
									 -c $(CORRESP) \
									 $${LOCAL:+-l}
	bash src/circos_conf/circos_input_gen.sh $(REF)


# Assembling transcripts and measuring coverage along genome
$(RNA)/assembled/ :
	bash $(RNA-SRC) -a $(BAM) -r $(OLD-REF) -o $@

# Contig correspondance file between old (unanchored)
# and new (anchored) assembly
$(CORRESP):
	bash src/convert_coord/corresp_contigs.sh -O $(OLD-REF) \
																						-N $(REF) \
																						-c $(CORRESP) \
																						2> $(LOG)/corresp.log

###############################################
#### 5. PI DIVERSITY FROM WILD SAMPLES WGS ####
###############################################

.PHONY : wgs_wild
wgs_wild : $(CORRESP) $(SIZES)
	#bash src/wgs_wild/bwa_wgs.sh --workdir $(WGS) --ref $(REF)
	#bash src/wgs_wild/qc_gen.sh --workdir $(WGS) --ref $(REF) --out $(WGS)/qc_output
	#bash src/wgs_wild/snps_wgs.sh --workdir $(WGS) --ref $(REF) --winsize 100
	Rscript src/wgs_wild/compute_PI.R -i $(WGS)/variant/chr.wild.matrix.txt \
									  -o $(WGS)/stats/win_w100_t10_PI.tsv \
									  -m 'window' --step_size 10 --win_size 100
	Rscript src/wgs_wild/compute_PI.R -i $(WGS)/variant/chr.wild.matrix.txt \
									  -o $(WGS)/stats/sites_PI.tsv -m 'site'
	bash src/convert_coord/CSD_contig.sh $(HITS) \
																			 3 \
																			 $(CORRESP) \
																			 $(ASSOC)

####################
#### MISC RULES ####
####################

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

# Check genome completeness using BUSCO
.PHONY : busco
busco :
	bash src/misc/genome_completeness.sh -r $(REF) \
																			 -d data/ref_genome/busco/


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
