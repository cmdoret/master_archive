include config.mk

.PHONY : all
all : $(MAP)/$(ALG)_$(MM)_$(K)_$(W)

# Running alignment with BWA
$(MAP)/$(ALG)_$(MM)_$(K)_$(W) : $(PROC)
	mkdir -p $@
	bsub bash src/mapping/bwa_script.sh $(ALG) $(MM) $(K) $(W)

# Running pstacks
#$(STACK) : $(MAP)/$(ALG)_$(MM)_$(K)_$(W)
#	bash src/stacks_pipeline/multi_pstacks.sh $<
#	for i in $(PROC)/*.fq.gz
#	do
#		bsub bash
#	done


.PHONY : clean
clean :
	rm -f demulti*
