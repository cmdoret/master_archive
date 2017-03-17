include config.mk

MM=2

$(wildcard $(PROC)/lib7/*) : $(wildcard $(RAW)/lib7/*)
	bash process_reads/sub_process.sh 7 $(MM)

$(wildcard $(PROC)/lib7b/*) : $(wildcard $(RAW)/lib7b/*)
	bash process_reads/sub_process.sh 7b $(MM)
