make_RL_combined_heatmap = function(dir,out_dir,threshold = 2.5, cell_type = NA )
{
  dir.create(out_dir,showWarnings = FALSE)
  
  ###################################
  # get name of all files to import #
  ###################################
  
  interaction_files = list.files(dir,"*interaction_scores.csv")
  receptor_files = list.files(dir,"*receptor_expression.csv")
  ligand_files = list.files(dir,"*ligand_expression.csv")
  
  
  num_files = length(interaction_files)
  
  
  ###########################
  # Read interaction scores #
  ###########################
  
  pair_names = list()
  cell_type_names = list()
  for (f in 1:num_files)
  {
    fName = ifelse(num_files==1, rev(unlist(strsplit(dir,"/",fixed=TRUE)))[1] # get model name from directory
                   , unlist(strsplit(interaction_files[f],"_",fixed = TRUE))[1]) # get model name file
    print(fName)
    
    interaction_scores = read.csv(paste(dir,interaction_files[f],sep=""),header=TRUE,row.names="Row")
    receptor_expression = read.csv(paste(dir,receptor_files[f],sep=""),header=TRUE)
    ligand_expression = read.csv(paste(dir,receptor_files[f],sep=""),header=TRUE)
    
    pair_names = union( pair_names, colnames( interaction_scores) )
    cell_type_names = union( cell_type_names, colnames(receptor_expression))
    
  }
  pair_names = sort(unlist(pair_names))
  cell_type_names = sort(unlist(cell_type_names))
  
  ####################################
  # add missing scores of cell-types #
  ####################################
  for (f in 1:num_files)
  {
    
    interaction_scores = read.csv(paste(dir,interaction_files[f],sep=""),header=TRUE,row.names="Row")
    receptor_expression = read.csv(paste(dir,receptor_files[f],sep=""),header=TRUE)
    ligand_expression = read.csv(paste(dir,receptor_files[f],sep=""),header=TRUE)
    
    # Interactions
    missing_interactions = setdiff(pair_names, colnames(interaction_scores) ) # Find interactions missing from model
    interaction_scores[missing_interactions] = 0 # set missing interaction scores to 0
    interaction_scores = interaction_scores[pair_names] # reorder columns 
    
    # Receptor
    missing_ct = setdiff(cell_type_names, colnames(receptor_expression) ) # Find interactions missing from model
    receptor_expression[missing_ct] = 0 # set missing interaction scores to 0
    receptor_expression = receptor_expression[cell_type_names] # reorder columns 
    # Ligand
    ligand_expression[missing_ct] = 0 # set missing interaction scores to 0
    ligand_expression = ligand_expression[cell_type_names] # reorder columns 
    
    if (f==1)
    { 
      combined_scores = interaction_scores 
      combined_receptor = receptor_expression
      combined_ligand = ligand_expression
    }
    else
    { 
      combined_scores = abind(combined_scores, interaction_scores,along = 3 )
      combined_receptor = abind(combined_receptor, receptor_expression,along = 3 )
      combined_ligand = abind(combined_ligand, ligand_expression,along = 3 )
    }
  }
  ###############################
  # Average across tumor models #
  ###############################
  interaction_scores = apply( combined_scores, MARGIN = c(1,2), FUN = mean)
  receptor_expression = apply( combined_receptor, MARGIN = c(1,2), FUN = mean)
  ligand_expression = apply( combined_ligand, MARGIN = c(1,2), FUN = mean)
  
  ####################################
  # Subset based on input cell-types #
  ####################################
  if (!any(is.na(cell_type)))
  {
    grep( paste(cell_type,collapse="|"),colnames(interaction_scores) )
    interaction_scores = interaction_scores[,grep( paste(cell_type,collapse="|"),colnames(interaction_scores) )]
    combined_scores = combined_scores[,grep( paste(cell_type,collapse="|"),colnames(combined_scores) ),]
  }
  
  ###################
  # Change row names#
  ###################
  split_names = strsplit(rownames(interaction_scores),"_") # Get receptor + ligand names from interaction names and find duplicates
  ligand_names = sapply( split_names, `[`,1)
  receptor_names = sapply( split_names, `[`,2)
  
  rownames(interaction_scores) = gsub("_","-", rownames(interaction_scores))
  colnames(interaction_scores) = gsub("_","-", colnames(interaction_scores))
  rownames(receptor_expression) = receptor_names
  rownames(ligand_expression) = ligand_names

  #########################################
  # Create receptor/ligand heatmap figure #
  #########################################
  
  # Display only rows/columns with any score greater than threshold
  row_idx = apply(interaction_scores>threshold,1,any)
  col_idx = apply(interaction_scores>threshold,2,any)

  receptor_subset = receptor_expression[row_idx,]
  ligand_subset = ligand_expression[row_idx,]


    
  ##### Set parameters for heatmaps
  col_func = colorRampPalette(c("white","red")) # Use red-white color scale for all heatmaps
  ht_global_opt(heatmap_column_names_gp = gpar(fontsize=12)) # font specs for heatmaps


  # Receptor expression
  duplicate_receptor = duplicated( receptor_subset ) 
  unique_receptor = receptor_subset[!duplicate_receptor,]
  rownames(unique_receptor) = rownames(receptor_subset)[!duplicate_receptor]
  
  cVals = seq( from = 0, to = ceiling(max(unique_receptor)), by = 0.02 )
  color_map = colorRamp2(cVals,col_func(length(cVals)))
  
  pdf(paste(out_dir,"sup_fig_3",cell_type,"_receptor.pdf"),width=11,height=8.5,paper='special') 
  hMap = Heatmap(unique_receptor, name = "receptor", column_title = "receptor",gap = unit(.02,"npc"),
                 row_names_side = "left",
                 na_col = "white", col = color_map, rect_gp = gpar(col = "gray"),
                 cluster_columns = FALSE,cluster_rows = FALSE,
                 show_heatmap_legend = TRUE)
  print(hMap)
  dev.off() # close pdf
  
  # Ligand expression
  duplicate_ligand = duplicated( rownames(ligand_subset) )
  unique_ligand = ligand_subset[!duplicate_ligand,]
  rownames(unique_ligand) = rownames(ligand_subset)[!duplicate_ligand]
  pdf(paste(out_dir,"sup_fig_3",cell_type,"_ligand.pdf"),width=11,height=8.5,paper='special') 
  
  cVals = seq( from = 0, to = ceiling(max(unique_ligand)), by = 0.02 )
  color_map = colorRamp2(cVals,col_func(length(cVals)))
  
  hMap = Heatmap(unique_ligand, name = "ligand", column_title = "ligand",gap = unit(.02,"npc"),
                 row_names_side = "left",
                 na_col = "white", col = color_map, rect_gp = gpar(col = "gray"),
                 cluster_columns = FALSE,cluster_rows = FALSE,
                 show_heatmap_legend = TRUE)
  print(hMap)
  dev.off() # close pdf

}