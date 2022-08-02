.PHONY: clean
.SHELL: /bin/bash

objects :=\
	data/processed/hic/YAPP_hic_loopCounts.rda\
	data/processed/hic/YAPP_hic_diff_loopCounts.rda\
	data/processed/microc/YAPP_microc_loopCounts.rda\
	data/processed/microc/YAPP_microc_diff_loopCounts.rda\
	plots/YAPP_diffLoops_sorb_PCA_hic.pdf\
	plots/YAPP_diffLoops_sorb_hcluster_hic.pdf\
	plots/YAPP_diffLoops_sorb_MA_hic.pdf\
	plots/YAPP_diffLoops_sorb_PCA_microc.pdf\
	plots/YAPP_diffLoops_sorb_hcluster_microc.pdf\
	plots/YAPP_diffLoops_sorb_MA_microc.pdf\

all: $(objects)

clean:
	rm -rf $(objects)

data/processed/hic/YAPP_hic_loopCounts.rda:\
	scripts/processing/genLoops_hic.R
	data/raw/hic/hg38/220722_dietJuicerCore/*
	data/raw/hic/hg38/sip-loops/isDroso/*/5kbLoops.txt
	data/raw/hic/hg38/sip-loops/noDroso/*/5kbLoops.txt
		mkdir -p data/processed/hic/
		Rscript scripts/processing/genLoops_hic.R
		
data/processed/hic/YAPP_hic_diff_loopCounts.rda:\
	scripts/processing/diffLoops_hic.R
	data/processed/hic/YAPP_hic_loopCounts.rda
		mkdir -p data/processed/hic/
		Rscript scripts/processing/diffLoops_hic.R

data/processed/microc/YAPP_microc_loopCounts.rda:\
	scripts/processing/genLoops_microc.R
	data/raw/microc/hg38_220801/220716_dietJuicerCore/*
	data/raw/microc/hg38_220801/sip-loops/*/5kbLoops.txt
		mkdir -p data/processed/microc/
		Rscript scripts/processing/genLoops_microc.R
		
data/processed/microc/YAPP_diff_loopCounts.rda:\
	scripts/processing/diffLoops_microc.R
	data/processed/microc/YAPP_microc_loopCounts.rda
		mkdir -p data/processed/microc/
		Rscript scripts/processing/diffLoops_microc.R
		
plots/YAPP_diffLoops_sorb_PCA_hic.pdf:\
	data/processed/hic/YAPP_hic_loopCounts.rda\
	scripts/processing/genLoops_hic.R
		mkdir -p plots
		Rscript scripts/genLoops_hic.R
	
plots/YAPP_diffLoops_sorb_hcluster_hic.pdf:\
	data/processed/hic/YAPP_hic_loopCounts.rda
	scripts/processing/genLoops_hic.R
	mkdir -p plots
	Rscript scripts/genLoops_hic.R

plots/YAPP_diffLoops_sorb_MA_hic.pdf:\
	data/processed/hic/YAPP_hic_loopCounts.rda
	scripts/processing/genLoops_hic.R
	mkdir -p plots
	Rscript scripts/genLoops_hic.R

plots/YAPP_diffLoops_sorb_PCA_microc.pdf:\
	data/processed/microc/YAPP_microc_loopCounts.rda
	scripts/processing/genLoops_microc.R
	mkdir -p plots
	Rscript scripts/genLoops_microc.R


plots/YAPP_diffLoops_sorb_hcluster_microc.pdf:\
	data/processed/microc/YAPP_microc_loopCounts.rda
	scripts/processing/genLoops_microc.R
	mkdir -p plots
	Rscript scripts/genLoops_microc.R

plots/YAPP_diffLoops_sorb_MA_microc.pdf:\
	data/processed/microc/YAPP_microc_loopCounts.rda
	scripts/processing/genLoops_microc.R
	mkdir -p plots
	Rscript scripts/genLoops_microc.R



		
	
	
	
	
	
	
	
	
	
	
		