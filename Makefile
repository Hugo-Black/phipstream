.PHONY: check pull test clean help

PYTHON ?= python3

help:
	@echo "Targets:"
	@echo "  check  verify the container engine and the Python dependency"
	@echo "  pull   fetch the pinned container images"
	@echo "  test   generate a synthetic dataset and run the pipeline over it"
	@echo "         set BEER=1 to include the enrichment step"
	@echo "  clean  remove results, work and log directories"

check:
	@bash bin/phipstream-setup

pull:
	@$(PYTHON) bin/containers.py pull

test: check
	@bash tests/end_to_end.sh

clean:
	rm -rf results/ work/ work_*/ logs/ .nextflow*
