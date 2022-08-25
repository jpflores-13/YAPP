.PHONY: clean
.SHELL: /bin/bash

objects :=\
	data/processed/hic/YAPP_hic_loopCounts.rds\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	data/processed/microc/YAPP_microc_loopCounts.rds\
	data/processed/microc/YAPP_microc_diff_loopCounts.rds\
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
	plots/meta_gainedLoops_APA.pdf\
	plots/meta_lostLoops_APA.pdf\
	plots/YAPP_HEK_hic_APA.pdf\
	plots/YAPP_HEK_microc_APA.pdf\
	external/HYPE/data/processed/hic/HYPE_loopCounts.rds\
	external/HYPE/data/processed/hic/HYPE_diff_loopCounts.rds\
	external/HYPE/plots/HYPE_diffLoops_nacl_PCA_hic.pdf\
	external/HYPE/plots/HYPE_diffLoops_nacl_hcluster_hic.pdf\
	external/HYPE/plots/HYPE_diffLoops_nacl_MA_hic.pdf\
	external/HYPE/plots/HYPE_T47D_hic_APA.pdf\
	vignettes/assets/meta_lostLoops_APA.png\
	vignetttes/assets/meta_gainedLoops_APA.png\
	vignettes/assets/HYPE_T47D_hic_APA.png\
	vignettes/assets/YAPP_HEK_microc_APA.png\
	vignettes/assets/YAPP_HEK_hic_APA.png\
	plots/known_motifEnrichment.pdf\
	plots/denovo_motifEnrichment.pdf
	

all: $(objects)
	echo done!

clean:
	rm -rf $(objects)
	
data/processed/hic/YAPP_hic_loopCounts.rds:\
	scripts/processing/genLoops_hic.R\
	data/raw/hic/hg38/220722_dietJuicerCore/*\
	data/raw/hic/hg38/sip-loops/isDroso/*/5kbLoops.txt\
	data/raw/hic/hg38/sip-loops/noDroso/*/5kbLoops.txt
		mkdir -p data/processed/hic/
		Rscript scripts/processing/genLoops_hic.R
		
data/processed/hic/YAPP_hic_diff_loopCounts.rds:\
	scripts/processing/diffLoops_hic.R\
	data/processed/hic/YAPP_hic_loopCounts.rds
		mkdir -p data/processed/hic/
		Rscript scripts/processing/diffLoops_hic.R
	
data/processed/microc/YAPP_microc_loopCounts.rds:\
	scripts/processing/genLoops_microc.R\
	data/raw/microc/hg38_220801/220716_dietJuicerCore/*\
	data/raw/microc/hg38_220801/sip-loops/isDroso/*/5kbLoops.txt\
	data/raw/microc/hg38_220801/sip-loops/noDroso/*/5kbLoops.txt
		mkdir -p data/processed/microc/
		Rscript scripts/processing/genLoops_microc.R
		
data/processed/microc/YAPP_microc_diff_loopCounts.rds:\
	scripts/processing/diffLoops_microc.R\
	data/processed/microc/YAPP_microc_loopCounts.rds
		mkdir -p data/processed/microc/
		Rscript scripts/processing/diffLoops_microc.R
		
plots/YAPP_diffLoops_sorb_PCA_hic.pdf:\
	data/processed/hic/YAPP_hic_loopCounts.rds\
	scripts/processing/genLoops_hic.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_hic.R
	
plots/YAPP_diffLoops_sorb_hcluster_hic.pdf:\
	data/processed/hic/YAPP_hic_loopCounts.rds\
	scripts/processing/genLoops_hic.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_hic.R

plots/YAPP_diffLoops_sorb_MA_hic.pdf:\
	data/processed/hic/YAPP_hic_loopCounts.rds\
	scripts/processing/genLoops_hic.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_hic.R

plots/YAPP_diffLoops_sorb_PCA_microc.pdf:\
	data/processed/microc/YAPP_microc_loopCounts.rds\
	scripts/processing/genLoops_microc.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_microc.R

plots/YAPP_diffLoops_sorb_hcluster_microc.pdf:\
	data/processed/microc/YAPP_microc_loopCounts.rds\
	scripts/processing/genLoops_microc.R
		mkdir -p plots
		Rscript scripts/processing/genLoops_microc.R

plots/YAPP_diffLoops_sorb_MA_microc.pdf:\
	data/processed/microc/YAPP_microc_loopCounts.rds\
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
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
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
	data/processed/microc/YAPP_microc_diff_loopCounts.rds\
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
	data/processed/microc/YAPP_microc_diff_loopCounts.rds\
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
	data/processed/microc/YAPP_microc_diff_loopCounts.rds\
	data/raw/microc/hg38_220801/sip-loops/*\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_omega/YAPP_HEK_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic\
	data/processed/microc/omega_bothDroso_loops.rds\
	data/processed/microc/sorb_bothDroso_loops.rds\
	data/processed/microc/cont_bothDroso_loops.rds
		mkdir -p plots
		Rscript scripts/analysis/YAPP_HEK_microc_lostLoops_rect.R
		
plots/meta_gainedLoops_APA.pdf:\
	scripts/analysis/meta_gainedLoops_APA.R\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/cont/HYPE_T47D_None_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/nacl/HYPE_T47D_NaCl_inter_30.hic
		mkdir -p plots
		Rscript scripts/analysis/meta_gainedLoops_APA.R
		
plots/meta_lostLoops_APA.pdf:\
	scripts/analysis/meta_lostLoops_APA.R\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/cont/HYPE_T47D_None_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/nacl/HYPE_T47D_NaCl_inter_30.hic
		mkdir -p plots
		Rscript scripts/analysis/meta_lostLoops_APA.R
		
plots/YAPP_HEK_hic_APA.pdf:\
	scripts/analysis/YAPP_HEK_hic_APA.R\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic
		mkdir -p plots
		Rscript scripts/analysis/YAPP_HEK_hic_APA.R

plots/YAPP_HEK_microc_APA.pdf:\
	scripts/analysis/YAPP_HEK_microc_APA.R\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic
		mkdir -p plots
		Rscript scripts/analysis/YAPP_HEK_microc_APA.R

external/HYPE/data/processed/hic/HYPE_loopCounts.rds:\
	external/HYPE/scripts/processing/genLoops.R\
	external/HYPE/data/raw/hic/hg38/220717_dietJuicerCore/*\
	external/HYPE/data/raw/hic/hg38/sip-loops/isDroso/*/5kbLoops.txt\
	external/HYPE/data/raw/hic/hg38/sip-loops/noDroso/*/5kbLoops.txt
		mkdir -p external/HYPE/data/processed/hic/
		Rscript external/HYPE/scripts/processing/genLoops.R
		
external/HYPE/data/processed/hic/HYPE_diff_loopCounts.rds:\
	external/HYPE/scripts/processing/diffAnalysis.R\
	external/HYPE/data/processed/hic/HYPE_loopCounts.rds
		mkdir -p external/HYPE/data/processed/hic/
		Rscript external/HYPE/scripts/processing/diffAnalysis.R
		
external/HYPE/plots/HYPE_diffLoops_nacl_PCA_hic.pdf:\
	external/HYPE/data/processed/hic/HYPE_loopCounts.rds\
	external/HYPE/scripts/processing/diffAnalysis.R
		mkdir -p external/HYPE/plots
		Rscript external/HYPE/scripts/processing/diffAnalysis.R
	
external/HYPE/plots/HYPE_diffLoops_nacl_hcluster_hic.pdf:\
	external/HYPE/data/processed/hic/HYPE_loopCounts.rds\
	external/HYPE/scripts/processing/diffAnalysis.R
		mkdir -p external/HYPE/plots
		Rscript external/HYPE/scripts/processing/diffAnalysis.R

external/HYPE/plots/HYPE_diffLoops_nacl_MA_hic.pdf:\
	external/HYPE/data/processed/hic/HYPE_loopCounts.rds\
	external/HYPE/scripts/processing/diffAnalysis.R
		mkdir -p external/HYPE/plots
		Rscript external/HYPE/scripts/processing/diffAnalysis.R

external/HYPE/plots/HYPE_hic_SCALE_tile_chr1.pdf:\
	external/HYPE/scripts/analysis/tileChr.R\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_omega/HYPE_T47D_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/cont/HYPE_T47D_None_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/nacl/HYPE_T47D_NaCl_inter_30.hic
		mkdir -p external/HYPE/plots
		Rscript external/HYPE/scripts/analysis/tileChr.R
		
external/HYPE/plots/HYPE_T47D_hic_APA.pdf:\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	external/HYPE/scripts/analysis/HYPE_T47D_hic_APA.R\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/cont/HYPE_T47D_None_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/nacl/HYPE_T47D_NaCl_inter_30.hic
		mkdir -p external/HYPE/plots
		Rscript external/HYPE/scripts/analysis/HYPE_T47D_hic_APA.R
		
vignettes/assets/meta_lostLoops_APA.png:\
	scripts/analysis/meta_lostLoops_APA.R
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/cont/HYPE_T47D_None_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/nacl/HYPE_T47D_NaCl_inter_30.hic
		mkdir -p vignettes/assets
		Rscript scripts/analysis/meta_lostLoops_APA.R

vignetttes/assets/meta_gainedLoops_APA.png:\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	data/raw/hic/hg38/subsamples_YAPP_HYPE/YAPP_HEK_control_subsampled_inter_30.hic\
	data/raw/hic/hg38/subsamples_YAPP_HYPE/YAPP_HEK_sorbitol_subsampled_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/cont/HYPE_T47D_None_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/nacl/HYPE_T47D_NaCl_inter_30.hic\
	scripts/analysis/meta_gainedLoops_APA.R
		mkdir -p vignettes/assets
		Rscript scripts/analysis/meta_gainedLoops_APA.R

vignettes/assets/HYPE_T47D_hic_APA.png:\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/cont/HYPE_T47D_None_inter_30.hic\
	external/HYPE/data/raw/hic/hg38/220718_dietJuicerMerge_treatment/nacl/HYPE_T47D_NaCl_inter_30.hic\
	external/HYPE/scripts/analysis/HYPE_T47D_hic_APA.R
		mkdir -p vignettes/assets
		Rscript external/HYPE/scripts/analysis/HYPE_T47D_hic_APA.R

vignettes/assets/YAPP_HEK_microc_APA.png:\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/cont/YAPP_HEK_cont_inter_30.hic\
	data/raw/microc/hg38_220801/220717_dietJuicerMerge_condition/sorb/YAPP_HEK_sorb_inter_30.hic\
	scripts/analysis/YAPP_HEK_microc_APA.R
		mkdir -p vignettes/assets
		Rscript scripts/analysis/YAPP_HEK_microc_APA.R

vignettes/assets/YAPP_HEK_hic_APA.png:\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/cont/YAPP_HEK_control_inter_30.hic\
	data/raw/hic/hg38/220716_dietJuicerMerge_condition/sorb/YAPP_HEK_sorbitol_inter_30.hic\
	scripts/analysis/YAPP_HEK_hic_APA.R
		mkdir -p vignettes/assets
		Rscript scripts/analysis/YAPP_HEK_hic_APA.R

data/processed/rna/YAPP_HEK_rnaseq_anchors.rds:\
	data/raw/rna/output/YAPP_HEK_WT_1_RNApipeSamplesheet.txt\
	scripts/processing/YAPP_HEK_rnaseq_anchors.R
		mkdir -p data/processed
		Rscript scripts/processing/YAPP_HEK_rnaseq_anchors.R

data/processed/rna/gainedLoop_genes.txt:\
	data/processed/rna/YAPP_HEK_rnaseq_anchors.rds\
	data/processed/hic/YAPP_hic_diff_loopCounts.rds\
	scripts/analysis/YAPP_HEK_rnaseq_anchors_boxplots.R
		mkdir -p data/processed
		Rscript scripts/analysis/YAPP_HEK_rnaseq_anchors_boxplots.RA
		
plots/known_motifEnrichment.pdf:\

plots/denovo_motifEnrichment.pdf:\







