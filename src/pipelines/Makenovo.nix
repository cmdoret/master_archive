include config.mk

.PHONY : all
all : $(POP)

# Building stacks denovo with Ustacks
$(USTACK) : $(PROC)
	rm -rf $@
	mkdir -p $@
	bash $(BWA-SRC) --mm $(MM) \
	                --ref $(REF) \
									--reads $(PROC) \
									--out $(MAP)

# Building catalog with Cstacks
#$(CSTACK): $(PSTACK)
	#rm -fr $@;
	#mkdir -p $@;
	#bash $(C-SRC) --pst $< \
	              --cst $@ \
								--lm $(LM)

# Map with GNSAP

# Integrate alignments with STACKS script
