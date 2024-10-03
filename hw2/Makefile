#phony targets 
.PHONY: carbon golf 

#flags corresponding to comand line arguments for all .py files in directory 
CARBON_PY=./carbon.py
GOLF_PY=./golf.py
WIDTH ?=
THETA ?=

#runs carbon.py script with provided argument
carbon:
	$(CARBON_PY) --plot=$(WIDTH)
#runs golf.py script with provided argument 
golf:
	$(GOLF_PY) --plot=$(THETA)

