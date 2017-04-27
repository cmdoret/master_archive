include config.mk

.PHONY : all
all : $(POP)

# Running alignment with BWA
$(MAP) : $(PROC)
	mkdir -p $@
	sed -i "s^\(main_dir=\).*^\1$(MAIN)^g" $(BWA-SRC);
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
	sed -i "s^\(wd=\).*^\1$(MAIN)/data^g" $(C-SRC);
	sed -i "s/\(MM=\)[0-9]*/\1$(LM)/g" $(C-SRC);
	bsub -K < $(C-SRC)

# Running Sstacks
$(SSTACK) : $(CSTACK)
	bash $(S-SRC) $<

# Running populations
$(POP) : $(SSTACK)
	mkdir -p $(POP)
	rm -rf $(DAT)/logs/populations
	mkdir -p $(DAT)/logs/populations
	sed -i "s^\(od=\).*^\1$(POP)^g" $(POP-SRC);
	sed -i "s^\(R=\).*^\1$(R)/data^g" $(POP-SRC);
	bsub -K $(POP-SRC)

.PHONY : clean
clean :
	rm -f *STDERR*
	rm -f *STDOUT*
	rm -f demulti*
	rm -rf bam
	rm -rf bsub_scripts
