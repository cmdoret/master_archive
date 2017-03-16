include config.mk

$(wildcard $(PROC)/lib7/*) : $(wildcard $(RAW)/lib7/*)
	bsub <./process_reads/process_7.sh

$(wildcard $(PROC)/lib7b/*) : $(wildcard $(RAW)/lib7b/*)
	bsub <./process_reads/process_7b.sh
