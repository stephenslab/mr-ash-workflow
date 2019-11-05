# A Makefile for generating the DSC results. Simply type "make" in
# this directory to run the DSC. Type "make clean" to clean up all the
# DSC outputs.

# This specifies the maximum number of threads used by dsc.
CORES=28

# RULES
# -----
all: mr_ash/mr_ash.db

mr_ash/mr_ash.db : mr_ash.dsc functions/*.R modules/simulate/*.R \
                   modules/fit/*.R modules/predict/*.R modules/score/*.R
	dsc -c $(CORES) mr_ash.dsc 

clean: 
	dsc --clean all mr_ash.dsc
