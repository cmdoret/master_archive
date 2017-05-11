include config.mk

.PHONY : all
all : $(POP)

# Running alignment with BWA
$(MAP) : $(PROC)
	rm -rf $@
	mkdir -p $@
	sed -i "s^\(main_dir=\).*^\1$(MAIN)^g" $(BWA-SRC)
	sed -i "s/\(MM=\)[0-9]*/\1$(MM)/g" $(BWA-SRC)
	sed -i "s/\(ALG=\)[a-z]*/\1$(ALG)/g" $(BWA-SRC)
	sed -i "s/\(K=\)[0-9]*/\1$(K)/g" $(BWA-SRC)
	sed -i "s/\(W=\)[0-9]*/\1$(W)/g" $(BWA-SRC)
	bsub -K <./$(BWA-SRC)

# Running Pstacks
$(PSTACK) : $(MAP)
	bash $(P-SRC) $< $(M)
	

# Running Cstacks
$(CSTACK) : $(PSTACK)
	rm -fr $@;
	mkdir -p $@;
	sed -i "s^\(wd=\).*^\1$(MAIN)/data^g" $(C-SRC)
	sed -i "s/\(MM=\)[0-9]*/\1$(LM)/g" $(C-SRC)
	bsub -K < $(C-SRC)

# Running Sstacks
$(SSTACK) : $(CSTACK)
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

.PHONY : lab_book
lab_book : $(LAB) $(MISC)
	rm  -f $(LAB)/*.log $(LAB)/*.synctex* $(LAB)/*.aux $(LAB)/*.out
	python $(MISC)/map_param.py
	bash $(MISC)/parse_pstacks.sh
	bash $(MISC)/parse_cstacks.sh
	bash $(MISC)/parse_VCF.sh
	Rscript $(MISC)/plot_VCF.R $(MISC)
	texi2pdf -b $(LAB)/lab_book.tex -c
	mv lab_book.pdf $(LAB)

.PHONY : clean
clean :
	rm -f *STDERR*
	rm -f *STDOUT*
	rm -f demulti*
	rm -rf bam
	rm -rf bsub_scripts

#.PHONY : ploidy
#ploidy:
#	bash $(MISC)/parse_VCF.sh
#	python src/ploidy/haplo_males.py $(VCFSUM)
