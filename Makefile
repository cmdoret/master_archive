include config.mk
RR=$(wildcard $(RAW)/lib*/)
LIB=$(patsubst $(RAW)/lib%/, %, $(RR))
PR=$(patsubst $(RAW)/lib%/, $(PROC)_$(ADAP_MM)/lib%/, $(RR))

.PHONY : run
run : $(PR)


# Processing raw reads of all libraries
$(PROC)_$(ADAP_MM)/lib%/ : process_reads/process_%.sh
	sed -i "s/MM=[0-9]/MM=$(ADAP_MM)/g" $<
	bsub <./$<


$(RAW)/qc/lib%.html : $(RAW)/lib%/)
	cat $(RAW)/lib$*/*.fq* > lib$*.fq
	mkdir -p $(RAW)/qc/
	sed -i "s/lib=[a-zA-Z0-9]*/lib=$*/g" process_reads/qc.sh
	bsub <./process_reads/qc.sh
	rm $(RAW)/lib$*/lib$*.fq

$(PROC)_$(MM)/qc/lib%.html : $(PROC)_$(MM)/lib%/
	cat $(PROC)/lib$*/*.fq* > lib$*.fq
	mkdir -p $(PROC)/qc/
	sed -i "s/lib=[a-zA-Z0-9]*/lib=$*/g" process_reads/post_qc.sh
	sed -i "s/MM=[0-9]/MM=$(ADAP_MM)/g" process_reads/post_qc.sh
	bsub <./post_qc.sh
	rm $(PROC)/lib$*/lib$*.fq
