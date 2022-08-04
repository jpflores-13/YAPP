.PHONY: clean
.SHELL: /bin/bash

objects :=\
	data/processed/hic/YAPP_hic_loopCounts.rda\
	data/processed/hic/YAPP_hic_diff_loopCounts.rda\
	data/processed/microc/YAPP_microc_loopCounts.rda\
	data/processed/microc/YAPP_microc_diff_loopCounts.rda\
	data/processed/hic/omega_bothDroso_loops.rds\
	data/processed/hic/cont_bothDroso_loops.rds\
	data/processed/hic/sorb_bothDroso_loops.rds\
	data/processed/microc/omega_bothDroso_loops.rds\
	data/processed/microc/cont_bothDroso_loops.rds\
	data/processed/microc/sorb_bothDroso_loops.rds\
	plots/YAPP_diffLoops_sorb_PCA_hic.pdf\
	plots/YAPP_diffLoops_sorb_hcluster_hic.pdf\
	plots/YAPP_diffLoops_sorb_MA_hic.pdf\
	plots/YAPP_diffLoops_sorb_PCA_microc.pdf\
	plots/YAPP_diffLoops_sorb_hcluster_microc.pdf\
	plots/YAPP_diffLoops_sorb_MA_microc.pdf\
	plots/YAPP_HEK_hic_gainedLoops_rect.pdf\
	plots/YAPP_HEK_hic_lostLoops_rect.pdf\
	plots/YAPP_HEK_microc_gainedLoops_rect.pdf\
	plots/YAPP_HEK_microc_lostLoops_rect.pdf\
	plots/YAPP_HEK_hic_APA.pdf\
	plots/YAPP_HEK_microc_APA.pdf

all: $(objects)
	echo done!

clean:
	rm -rf $(objects)
	
data/processed/hic/YAPP_hic_loopCounts.rda:\
	scripts/processing/genLoops_hic.R\
	data/raw/hic/hg38/220722_dietJuicerCore/*\
	data/raw/hic/hg38/sip-loops/isDroso/*/5kbLoops.txt\
	data/raw/hic/hg38/sip-loops/noDroso/*/5kbLoops.txt
		mkdir -p data/processed/hic/
		Rscript scripts/processing/genLoops_hic.R
		
data/processed/hic/YAPP_hic_diff_loopCounts.rda:\
	scripts/processing/diffLoops_hic.R\
	data/processed/hic/YAPP_hic_loopCounts.rda
		mkdir -p data/processed/hic/
		Rscript scripts/processing/diffLoops_hic.R
		
data/processed/hic/omega_bothDroso_loops.rds:\
	data/raw/hic/hg38/sip-loops/isDroso/omega/5kbLoops.txt\
	data/raw/hic/hg38/sip-loops/noDroso/omega/5kbLoops.txt\
	scripts/processing/loopLists_hic_condition.R
		mkdir -p data/processed/hic/
		Rscript scripts/processing/loopLists_hic_condition.R
	
data/processed/hic/sorb_bothDroso_loops.rds:\
	data/raw/hic/hg38/sip-loops/isDroso/sorb/5kbLoops.txt\
	data/raw/hic/hg38/sip-loops/noDroso/sorb/5kbLoops.txt\
	scripts/processing/loopLists_hic_condition.R
		mkdir -p data/processed/hic/
		Rscript scripts/processing/loopLists_hic_condition.R
	
data/processed/hic/cont_bothDroso_loops.rds:\
	data/raw/hic/hg38/sip-loops/isDroso/cont/5kbLoops.txt\
	data/raw/hic/hg38/sip-loops/noDroso/cont/5kbLoops.txt\
	scripts/processing/loopLists_hic_condition.R
		mkdir -p data/processed/hic/
		Rscript scripts/processing/loopLists_hic_condition.R		
	
data/processed/microc/omega_bothDroso_loops.rds:\
	data/raw/microc/hg38_220801/sip-loops/isDroso/omega/5kbLoops.txt\
	data/raw/microc/hg38_220801/sip-loops/noDroso/omega/5kbLoops.txt\
	scripts/processing/loopLists_microc_condition.R
		mkdir -p data/processed/microc/
		Rscript scripts/processing/loopLists_microc_condition.R
	
data/processed/microc/sorb_bothDroso_loops.rds:\
	data/raw/microc/hg38_220801/sip-loops/isDroso/sorb/5kbLoops.txt\
	data/raw/microc/hg38_220801/sip-loops/noDroso/sorb/5kbLoops.txt\
	scripts/processing/loopLists_microc_condition.R
		mkdir -p data/processed/microc/
		Rscript scripts/processing/loopLists_microc_condition.R
	
data/processed/microc/cont_bothDroso_loops.rds:\
	scripts/processing/loopLists_microc_condition.R
		mkdir -p data/processed/microc/
		Rscript scripts/processing/loopLists_microc_condition.R		

data/processed/microc/YAPP_microc_loopCounts.rda:\
	scripts/processing/genLoops_microc.R\
	data/raw/microc/hg38_220801/220716_dietJuicerCore/*\
	data/raw/microc/hg38_220801/sip-loops/isDroso/*/5kbLoops.txt\
	data/raw/microc/hg38_220801/sip-loops/noDroso/*/5kbLoops.txt
		mkdir -p data/processed/microc/
		Rscript scripts/processing/genLoops_microc.R
		
data/processed/microc/YAPP_microc_diff_loopCounts.rda:\
	scripts/processing/diffLoops_microc.R\
	data/processed/microc/YAPP_microc_loopCounts.rda
		mkdir -p data/processed/microc/
		Rscript scripts/processing/diffLoops_microc.R
		
plots/YAPP_diffLoops_sorb_PCA_hic.pdf:\
	data/processed/hic/YAPP_hic_loopCounts.rda\
	scripts/processing/genLoops_hic.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_hic.R
	
plots/YAPP_diffLoops_sorb_hcluster_hic.pdf:\
	data/processed/hic/YAPP_hic_loopCounts.rda\
	scripts/processing/genLoops_hic.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_hic.R

plots/YAPP_diffLoops_sorb_MA_hic.pdf:\
	data/processed/hic/YAPP_hic_loopCounts.rda\
	scripts/processing/genLoops_hic.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_hic.R

plots/YAPP_diffLoops_sorb_PCA_microc.pdf:\
	data/processed/microc/YAPP_microc_loopCounts.rda\
	scripts/processing/genLoops_microc.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_microc.R

plots/YAPP_diffLoops_sorb_hcluster_microc.pdf:\
	data/processed/microc/YAPP_microc_loopCounts.rda\
	scripts/processing/genLoops_microc.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_microc.R

plots/YAPP_diffLoops_sorb_MA_microc.pdf:\
	data/processed/microc/YAPP_microc_loopCounts.rda\
	scripts/processing/genLoops_microc.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_microc.R

plots/YAPP_hic_SCALE_tile_chr1.pdf:\
	scripts/analysis/tileChr_hic.R
		mkdir -p plots
		Rscript scripts/analysis/tileChr_hic.R
		
plots/YAPP_microc_SCALE_tile_chr1.pdf:\
	scripts/analysis/tileChr_microc.R
		mkdir -p plots 
		Rscript scripts/analysis/tileChr_microc.R
		
plots/YAPP_HEK_hic_gainedLoops_rect.pdf:\
	scripts/analysis/YAPP_HEK_hic_gainedLoops_rect.R\
	data/processed/hic/YAPP_hic_diff_loopCounts.rda\
	data/raw/hic/hg38/sip-loops/*\
	data/raw/hic/hg38/220716_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic\
	data/processed/hic/omega_bothDroso_loops.rds\
	data/processed/hic/sorb_bothDroso_loops.rds\
	data/processed/hic/cont_bothDroso_loops.rds
		mkdir -p plots
		Rscript scripts/analysis/YAPP_HEK_hic_gainedLoops_rect.R
		
plots/YAPP_HEK_hic_lostLoops_rect.pdf:\
	scripts/analysis/YAPP_HEK_hic_lostLoops_rect.R\
	data/processed/microc/YAPP_microc_diff_loopCounts.rda\
	data/raw/hic/hg38/sip-loops/*\
	data/raw/hic/hg38/220716_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic\
	data/processed/hic/omega_bothDroso_loops.rds\
	data/processed/hic/sorb_bothDroso_loops.rds\
	data/processed/hic/cont_bothDroso_loops.rds
		mkdir -p plots
		Rscript scripts/analysis/YAPP_HEK_hic_lostLoops_rect.R
		
plots/YAPP_HEK_microc_gainedLoops_rect.pdf:\
	scripts/analysis/YAPP_HEK_microc_gainedLoops_rect.R\
	data/processed/microc/YAPP_microc_diff_loopCounts.rda\
	data/raw/microc/hg38_220801/sip-loops/*\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic\
	data/processed/microc/omega_bothDroso_loops.rds\
	data/processed/microc/sorb_bothDroso_loops.rds\
	data/processed/microc/cont_bothDroso_loops.rds
		mkdir -p plots
		Rscript scripts/analysis/YAPP_HEK_microc_gainedLoops_rect.R	
	
plots/YAPP_HEK_microc_lostLoops_rect.pdf:\
	scripts/analysis/YAPP_HEK_microc_lostLoops_rect.R\
	data/processed/microc/YAPP_microc_diff_loopCounts.rda\
	data/raw/microc/hg38_220801/sip-loops/*\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic\
	data/processed/microc/omega_bothDroso_loops.rds\
	data/processed/microc/sorb_bothDroso_loops.rds\
	data/processed/microc/cont_bothDroso_loops.rds
		mkdir -p plots
		Rscript scripts/analysis/YAPP_HEK_microc_lostLoops_rect.R
		
plots/YAPP_HEK_hic_APA.pdf:\
	scripts/analysis/YAPP_HEK_hic_APA.R\
	data/processed/hic/YAPP_hic_loopCounts.rda\
	data/processed/hic/omega_bothDroso_loops.rds\
	data/processed/hic/cont_bothDroso_loops.rds\
	data/processed/hic/sorb_bothDroso_loops.rds\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic
		mkdir -p plots
		Rscript scripts/analysis/YAPP_HEK_hic_APA.R
		
plots/YAPP_HEK_microc_APA.pdf:\
	scripts/analysis/YAPP_HEK_microc_APA.R\
	data/processed/microc/YAPP_microc_loopCounts.rda\
	data/processed/microc/omega_bothDroso_loops.rds\
	data/processed/microc/cont_bothDroso_loops.rds\
	data/processed/microc/sorb_bothDroso_loops.rds\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic
		mkdir -p plots
		Rscript scripts/analysis/YAPP_HEK_microc_APA.R

	
	
	
	
	
		