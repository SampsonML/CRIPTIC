# Makefile for criptic
.PHONY: all clean

PROB ?= NONE

all: criptic

debug: criptic-debug

criptic: setprob
	cd Source && $(MAKE)
	cp Source/criptic Run/

clean:
	cd Source && $(MAKE) clean

setprob:
ifneq ($(PROB),NONE)
	@echo "*** Setting problem to $(PROB) ***"
	mkdir Run
	cp -f Prob/$(PROB)/Prob.cpp Source/
	cp -f Prob/$(PROB)/Definitions.H Source/
	cp -f Prob/$(PROB)/criptic.in Run/
	cp -f Prob/$(PROB)/Input_01 Run/
endif
