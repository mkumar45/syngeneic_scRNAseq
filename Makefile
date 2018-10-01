# Makefile for cell-cell communication analysis

.PHONY : help
help : Makefile
	@sed -n 's/^##//p' $<

.PHONY : all
all : results/syngeneic_tsne.csv
	

# $@ = target of the current rule
# $^ = all dependencies
# $< = first dependices of current rule
# % = wildcard (use in dependencies)
# $* = replaced by the stem the rule matched (use in actions)


################################
### Syngeneic mouse analysis ###
################################							  
##
results/syngeneic_processed_expression.csv : data/syngeneic_expression.csv data/syngeneic_row_data.csv data/HK_genes.txt code/helper/find_homologs.m code/preprocessing/preprocess_syngeneic.m
	matlab -r "addpath(genpath('code/')); preprocess_syngeneic( 'data/syngeneic_expression.csv', 'data/syngeneic_row_data.csv', 'data/syngeneic_column_data.csv', \
																1000, '$@'); quit()"

## results/syngeneic_tsne.csv : Run t-SNE on processed syngeneic expression and save output to file
results/syngeneic_tsne.csv : code/run_tsne.m code/bhtsne/bh_tsne results/syngeneic_processed_expression.csv
	matlab -r "addpath(genpath('code/'));run_tsne( 'results/syngeneic_processed_expression.csv', '$@');quit()"

	  
## code/bhtsne/bh_tsne : Compile C++ bhtsne implementation into excucable for Linux. Must be done before running t-SNE. See bhtsne README.
code/bhtsne/bh_tsne : code/bh_tsne.tar.gz
	tar -xzf code/bh_tsne.tar.gz --directory code/
	g++ code/bhtsne/sptree.cpp code/bhtsne/tsne.cpp code/bhtsne/tsne_main.cpp -o code/bhtsne/bh_tsne

## results/%_classification/ : predict cell-types by training classifier using using training data created using pre-defined cell-type markers
results/syngeneic_classification/ :	code/cell_type_classification.m results/syngeneic_processed_expression.csv data/syngeneic_cell_type_markers.csv data/syngeneic_row_data.csv
	mkdir -p $@
	matlab -r "addpath(genpath('code/')); \
				cell_type_classification( 'results/syngeneic_processed_expression.csv' ,'data/syngeneic_row_data.csv', 'data/syngeneic_cell_type_markers.csv', '$@' ); \
				quit()"
	
results/syngeneic_classification/cell_type_predictions.csv : results/syngeneic_classification/

	
## data/mouse_receptor_ligand.csv : convert receptor-ligand pairs from Ramilowski et al. to mouse gene symobls	
data/mouse_receptor_ligand.csv : data/receptor_ligand_Ramilowski_modified.csv code/helper/make_mouse_RL.m
	matlab -r "addpath(genpath('code/')); make_mouse_RL( '$<'); quit()"

			
## results/syngeneic_communication/ : calculate scores for receptor-ligand interactions between all pairs of classified cell-types				
results/syngeneic_communication/ : code/cell_cell_communication.m code/communication/ \
						   results/syngeneic_processed_expression.csv data/syngeneic_row_data.csv results/syngeneic_processed_column_data.csv data/mouse_receptor_ligand.csv \
						   results/syngeneic_classification/cell_type_predictions.csv	
	mkdir -p $@
	matlab -r "addpath(genpath('code/')); cell_cell_communication('results/syngeneic_processed_expression.csv', \
																  'data/syngeneic_row_data.csv', \
																  'results/syngeneic_processed_column_data.csv', \
																  'results/syngeneic_classification/cell_type_predictions.csv', \
																  'data/mouse_receptor_ligand.csv', \
																  'Model', \
																  'results/syngeneic_communication/'); quit()"
																  

results/syngeneic_phenotype.csv : code/phenotype/make_phenotype_table.m results/syngeneic_processed_column_data.csv results/syngeneic_classification/cell_type_predictions.csv data/syngeneic_tumor_growth.xlsx
	matlab -r "addpath(genpath('code/')); make_phenotype_table( 'results/syngeneic_classification/cell_type_predictions.csv', \
																'results/syngeneic_processed_column_data.csv', \
																'data/syngeneic_tumor_growth.xlsx', \
																'$@' ); quit()"
																
results/syngeneic_phenotype/ : code/phenotype/phenotype_correlation.m results/syngeneic_phenotype.csv
	mkdir -p $@
	matlab -r "addpath(genpath('code/')); phenotype_correlation( 'results/syngeneic_communication/', \
																 'results/syngeneic_phenotype.csv', \
																  100, 'results/syngeneic_phenotype/' ); quit()"

																  
																  
#################################
### Human metastatic melanoma ###
#################################
																  
## results/tirosh_processed_expression.csv results/tirosh_row_data.csv results/tirosh_column_data.csv : preprocess raw data from Tirosh et al.
results/tirosh_processed_expression.csv results/tirosh_row_data.csv results/tirosh_column_data.csv : data/GSE72056_melanoma_single_cell_revised_v2.txt code/preprocessing/preprocess_tirosh.m
	matlab -r "addpath(genpath('code/')); preprocess_tirosh( '$<', 1500 ); quit()"																  
																  											  
## results/tirosh_tsne.csv : Run t-SNE on processed tirosh expression and save output to file
results/tirosh_tsne.csv : results/tirosh_processed_expression.csv code/run_tsne.m code/bhtsne/bh_tsne 
	matlab -r "addpath(genpath('code/')); run_tsne( '$<', '$@');quit()"
	
	
## results/%_classification/ : predict cell-types by training classifier using using training data created using pre-defined cell-type markers
results/tirosh_classification/ : results/tirosh_processed_expression.csv data/tirosh_cell_type_markers.csv results/tirosh_processed_row_data.csv code/cell_type_classification.m
	mkdir -p $@
	matlab -r "addpath(genpath('code/')); \
				cell_type_classification( '$<' ,'results/tirosh_processed_row_data.csv', 'data/tirosh_cell_type_markers.csv', '$@' ); \
				quit()"
				
## results/%_classification/ : predict cell-types by training classifier using using training data created using pre-defined cell-type markers
results/tirosh_classification/tcell_subsets/ : results/tirosh_processed_expression.csv results/tirosh_classification/cell_type_predictions.csv data/tirosh_tcell_markers.csv results/tirosh_processed_row_data.csv code/tirosh/tcell_classification.m
	mkdir -p $@
	matlab -r "addpath(genpath('code/')); \
				tcell_classification( '$<' , 'results/tirosh_classification/cell_type_predictions.csv','results/tirosh_processed_row_data.csv', 'data/tirosh_tcell_markers.csv', '$@', 0.95 ); \
				quit()"
												  
## results/tirosh_communication/ : calculate scores for receptor-ligand interactions between all pairs of classified cell-types				
results/tirosh_communication/ : code/tirosh/tirosh_communication.m code/communication/ \
						   results/tirosh_processed_expression.csv results/tirosh_processed_row_data.csv results/tirosh_processed_column_data.csv data/receptor_ligand_Ramilowski_modified.csv \
						   results/tirosh_classification/tcell_subsets/cell_type_labels.txt	
	mkdir -p $@
	matlab -r "addpath(genpath('code/')); tirosh_communication('results/tirosh_processed_expression.csv', \
																  'results/tirosh_processed_row_data.csv', \
																  'results/tirosh_processed_column_data.csv', \
																  'results/tirosh_classification/tcell_subsets/cell_type_labels.txt', \
																  'data/receptor_ligand_Ramilowski_modified.csv', \
																  'tumor_id', \
																  '$@'); quit()"																  
															  
## results/tirosh_phenotype.csv :						  
results/tirosh_phenotype.csv : results/tirosh_processed_expression.csv code/phenotype/make_tirosh_phenotype.m results/tirosh_processed_column_data.csv
	matlab -r "addpath(genpath('code/')); make_tirosh_phenotype( '$<', \
																 'results/tirosh_classification/tcell_subsets/cell_type_labels.txt', \
																 'results/tirosh_processed_column_data.csv', \
																 'results/tirosh_processed_row_data.csv', \
																 '$@' ); quit()"
																 
																 
																 ## results/tirosh_phenotype.csv :						  
results/tirosh_phenotype/ : code/phenotype/phenotype_correlation.m results/tirosh_phenotype.csv
	mkdir -p $@
	matlab -r "addpath(genpath('code/')); phenotype_correlation( 'results/tirosh_communication/', \
																 'results/tirosh_phenotype.csv', \
																  100, '$@' ); quit()"
###############
### Figures ###
###############															  
																  
																  
## figures : create main and/or supplemental figures
.PHONY : figures_all figure1
figures_all : figure1

figure1 : code/plotting/make_figure_1.m results/syngeneic_tsne.csv results/syngeneic_classification/cell_type_predictions.csv
	matlab -r "addpath(genpath('code/')); \
			   make_figure_1( 'results/syngeneic_tsne.csv', 'results/syngeneic_processed_column_data.csv', \
			   'results/syngeneic_classification/cell_type_predictions.csv', 'figures/' ); quit()"
	
figure2 : code/plotting/make_score_heatmap.R code/plotting/read_null_scores.R results/syngeneic_communication/
	Rscript code/plotting/make_score_heatmap.R results/syngeneic_communication/ figures/ 1
## clean : remove processed files
.PHONY : clean_all clean_processed clean_tsne clean_classification clean_communicatoin
clean_all : clean_processed clean_tsne clean_classification clean_communication

clean_processed : 
	rm data/*_processed.csv
clean_tsne :
	rm results/*_tsne.csv
clean_classification :
	rm results/*_classification/