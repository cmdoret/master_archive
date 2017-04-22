include config.mk

.PHONY : all
all : $(CSTACK)

# Running alignment with BWA
$(MAP) : $(PROC)
	mkdir -p $@
	sed -i "s^\(main_dir=\).*^\1$(MAIN)^g" $(BWA-SRC);
	sed -i "s/\(MM=\)[0-9]*/\1$(MM)/g" $(BWA-SRC)
	sed -i "s/\(ALG=\)[a-z]*/\1$(ALG)/g" $(BWA-SRC)
	sed -i "s/\(K=\)[0-9]*/\1$(K)/g" $(BWA-SRC)
	sed -i "s/\(W=\)[0-9]*/\1$(W)/g" $(BWA-SRC)
	bsub -K <./$(BWA-SRC)

# Running pstacks
$(PSTACK) : $(MAP)
	bash $(P-SRC) $< $(M)

# Running cstacks
$(CSTACK) : $(PSTACK)
	rm -fr $@;
	mkdir -p $@;
	sed -i "s^\(wd=\).*^\1$(MAIN)/data^g" $(C-SRC);
	sed -i "s/\(MM=\)[0-9]*/\1$(LM)/g" $(C-SRC);
	bsub -K < $(C-SRC)

.PHONY : clean
clean :
	rm -f *STDERR*
	rm -f *STDOUT*
	rm -f demulti*
	rmdir bam
