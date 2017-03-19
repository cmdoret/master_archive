include config.mk
RR=$(wildcard $(RAW)/lib*)
LIB=$(patsubst $(RAW)/lib%/, %, $(RR))
PR=$(patsubst $(RAW)/lib%/, $(PROC)_$(ADAP_MM)/lib%/, $(RR))
PQC_FILES=$(patsubst $(RAW)/lib%, $(PROC)_$(ADAP_MM)/qc/lib%.html, $(RR))

.PHONY : all
all : $(PQC_FILES)

# Processing raw reads of all libraries
$(PROC)_$(ADAP_MM)/lib%/ : src/process_reads/process_%.sh
	mkdir -p $@
	sed -i "s/MM=[0-9]/MM=$(ADAP_MM)/g" $<
	bsub <./$<


$(RAW)/qc/lib%.html : $(RAW)/lib%/
	rm $</lib$*.fq
	cat $</*.fq* > $</lib$*.fq
	mkdir -p $(RAW)/qc/
	sed -i "s/lib=[a-zA-Z0-9]*/lib=$*/g" src/process_reads/qc.sh
	bsub <./src/process_reads/qc.sh

$(PROC)_$(ADAP_MM)/qc/lib%.html : $(PROC)_$(ADAP_MM)/lib%/
	rm -f $</lib$*.fq
	cat $</*.fq* > $</lib$*.fq
	mkdir -p $(PROC)_$(ADAP_MM)/qc/
	sed -i "s/lib=[a-zA-Z0-9]*/lib=$*/g" src/process_reads/post_qc.sh
	sed -i "s/MM=[0-9]/MM=$(ADAP_MM)/g" src/process_reads/post_qc.sh
	bsub <./src/process_reads/post_qc.sh

.PHONY : clean
clean : 
	rm -f demulti*
