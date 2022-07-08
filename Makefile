.PHONY: clean
.SHELL: /bin/bash

objects :=\
	data/processed/hic/isDroso/loopCounts.rda\
	data/processed/hic/isDroso/diff_loopCounts.rda\
	data/processed/hic/noDroso/loopCounts.rda\
	data/processed/hic/noDroso/diff_loopCounts.rda\
	plots/YAPP_diffLoops_sorb_MA_isDroso_betaPriorT.pdf

all: $(objects)

clean:
	rm -rf $(objects)

data/processed/hic/isDroso/loopCounts.rda:\
	scripts/processing/isDroso/genLoops_isDroso.R
	data/raw/hic/hg38/220627_dietJuicerCore/*
	data/raw/hic/hg38/sip-loops/isDroso/*/5kbLoops.txt
		mkdir -p data/processed/hic/isDroso
		Rscript scripts/processing/isDroso/genLoops_isDroso.R
		
data/processed/hic/isDroso/diff_loopCounts.rda:\
	scripts/processing/isDroso/genLoops_isDroso.R
	data/processed/hic/isDroso/loopCounts.rda
		mkdir -p data/processed/hic/isDroso/
		Rscript scripts/processing/isDroso/diffAnalysisYAPP_isDroso.R

data/processed/hic/noDroso/loopCounts.rda:\
	scripts/processing/noDroso/genLoops_noDroso.R
	data/raw/hic/hg38/220627_dietJuicerCore/*
	data/raw/hic/hg38/sip-loops/noDroso/*/5kbLoops.txt
		mkdir -p data/processed/hic/noDroso
		Rscript scripts/processing/noDroso/genLoops_noDroso.R
		
data/processed/hic/noDroso/diff_loopCounts.rda:\
	scripts/processing/noDroso/genLoops_noDroso.R
	data/processed/hic/noDroso/loopCounts.rda
		mkdir -p data/processed/hic/noDroso
		Rscript scripts/processing/noDroso/diffAnalysisYAPP_noDroso.R
		
plots/YAPP_diffLoops_sorb_MA_isDroso_betaPriorT.pdf:\
	scripts/processing/isDroso/diffAnalysisYAPP_isDroso.R
	data/processed/hic/isDroso/loopCounts.rda
		mkdir -p plots
		Rscript scripts/processing/isDroso/diffAnalysisYAPP_isDroso.R
	
	
	
	
	
	
	
	
	
	
		