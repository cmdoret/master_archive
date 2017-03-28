include config.mk

.PHONY : all
all : $(MAP)/$(ALG)_$(MM)_$(K)_$(W)

# Running alignment with BWA
$(MAP)/$(ALG)_$(MM)_$(K)_$(W) : $(PROC)
	mkdir -p $@
	sed -i "s/\(MM=\)[0-9]*/\1$(MM)/g" src/mapping/bwa_script.sh
	sed -i "s/\(ALG=\)[a-z]*/\1$(ALG)/g" src/mapping/bwa_script.sh
	sed -i "s/\(K=\)[0-9]*/\1$(K)/g" src/mapping/bwa_script.sh
	sed -i "s/\(W=\)[0-9]*/\1$(W)/g" src/mapping/bwa_script.sh
	bsub <./src/mapping/bwa_script.sh

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
	rmdir bam
