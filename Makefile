include config.mk

# Processing reads of library 7
$(wildcard $(PROC)/lib7/*) : $(wildcard $(RAW)/lib7/*) process_reads/process_7.sh
	#bash process_reads/sub_process.sh 7 $(ADAP_MM)
	sed -i "s/MM=[0-9]/MM=$(ADAP_MM)/g" process_reads/process_7.sh
	bsub <./process_reads/process_7.sh

# Processing reads of library 7b
$(wildcard $(PROC)/lib7b/*) : $(wildcard $(RAW)/lib7b/*) process_reads/process_7b.sh
	#bash process_reads/sub_process.sh 7b $(ADAP_MM)
	sed -i "s/MM=[0-9]/MM=$(ADAP_MM)/g" process_reads/process_7b.sh
	bsub <./process_reads/process_7b.sh

.PHONY : qc
qc: $(wilcard $(RAW)/lib%/*)
	module add UHTS/Quality_control/fastqc/0.11.2
	cat $(RAW)/lib$*/*.fq* > lib$*.fq
	mkdir -p $(RAW)/qc/
	bsub -q normal fastqc -o $(RAW)/qc/ $(RAW)/lib$*/lib$*.fq
	module rm UHTS/Quality_control/fastqc/0.11.2
	rm $(RAW)/lib$*/lib$*.fq

.PHONY : post_qc
post_qc: $(wilcard $(PROC)/lib%/*)
	module add UHTS/Quality_control/fastqc/0.11.2
	cat $(PROC)/lib$*/*.fq* > lib$*.fq
	mkdir -p $(PROC)/qc/
	bsub -q normal fastqc -o $(PROC)/qc/ $(PROC)/lib$*/lib$*.fq
	module rm UHTS/Quality_control/fastqc/0.11.2
	rm $(PROC)/lib$*/lib$*.fq
