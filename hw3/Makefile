#phony targets 
.PHONY: oscillator poisson

#flags corresponding to command line input
OS_PY=./oscillator.py 
PO_PY=./poisson.py
PART ?= 

#runs oscillator script with arguments:
osci:
	$(OS_PY) --part=$(PART)
#runs poisson script with arguments:
pois:
	$(PO_PY) --part=$(PART)
