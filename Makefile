DIRTY_FILES:=$(shell find . -name "*.pyc") $(shell find . -name "*.pickle")
DIRTY_FILES:=$(DIRTY_FILES) $(shell find . -name "*.aux")
DIRTY_FILES:=$(DIRTY_FILES) $(shell find . -name "*.toc")
DIRTY_FILES:=$(DIRTY_FILES) $(shell find . -name "*.log")
DIRTY_FILES:=$(DIRTY_FILES) $(shell find . -name "*.pdf")

DIRTY_DIRS=$(shell find . -name "__pycache__")

.PHONY: clean
clean:
	rm -rf $(DIRTY_DIRS)
	rm -f $(DIRTY_FILES)

.PHONY: doc
doc:
	cd doc; latex -output-format=pdf notes.tex
