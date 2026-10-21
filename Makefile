PYTHON ?= python3
CONFIG ?= configs/full.json
OUT ?= results/$(notdir $(basename $(CONFIG)))
export MPLCONFIGDIR := $(CURDIR)/.cache/matplotlib

.PHONY: all simulate analyze figures clean test

all: figures

simulate:
	$(PYTHON) -m src.simulation --config "$(CONFIG)" --output "$(OUT)"

analyze: simulate
	$(PYTHON) -m src.analyze --config "$(CONFIG)" --output "$(OUT)"

figures: analyze
	$(PYTHON) -m src.figures --config "$(CONFIG)" --output "$(OUT)"

test:
	$(PYTHON) -m pytest -q

clean:
	$(PYTHON) -m src.clean --config "$(CONFIG)" --output "$(OUT)"
