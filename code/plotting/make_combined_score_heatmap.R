make_combined_score_heatmap = function(dir,out_dir,threshold = 2.5, cell_type = NA, FDR = 0.1, normalize = identity )
{
  # 
  file_list = list.files(dir,"*interaction_scores.csv")
  num_files = length(file_list)
  dir.create(out_dir,showWarnings = FALSE)
  
  
  
  ##### Set parameters for heatmaps
  col_func = colorRampPalette(c("white","red")) # Use red-white color scale for all heatmaps
  ht_global_opt(heatmap_column_names_gp = gpar(fontsize=12)) # font specs for heatmaps
  
  ### Read interaction scores
  interaction_names = list()
  for (f in 1:length(file_list))
  {
    fName = ifelse(num_files==1, rev(unlist(strsplit(dir,"/",fixed=TRUE)))[1] # get model name from directory
                   , unlist(strsplit(file_list[f],"_",fixed = TRUE))[1]) # get model name file
    print(fName)
    interaction_scores = read.csv(paste(dir,file_list[f],sep=""),header=TRUE,row.names="Row")
    interaction_names = union( interaction_names, colnames( interaction_scores) )
  }
  interaction_names = sort(unlist(interaction_names))
  
  for (f in 1:length(file_list))
  {
    fName = ifelse(num_files==1, rev(unlist(strsplit(dir,"/",fixed=TRUE)))[1] # get model name from directory
                   , unlist(strsplit(file_list[f],"_",fixed = TRUE))[1]) # get model name file
    
    interaction_scores = read.csv(paste(dir,file_list[f],sep=""),header=TRUE,row.names="Row")
    
    missing_interactions = setdiff(interaction_names, colnames(interaction_scores) ) # Find interactions missing from model
    interaction_scores[missing_interactions] = 0 # set missing interaction scores to 0
    interaction_scores = interaction_scores[interaction_names] # reorder columns 
    
    if (f==1)
      { combined_scores = interaction_scores }
    else
      { combined_scores = abind(combined_scores, interaction_scores,along = 3 )}
  }
  
  interaction_scores = apply( combined_scores, MARGIN = c(1,2), FUN = mean)
  
  ##### Subset columns based on input cell-type to display #
  if (!any(is.na(cell_type)))
  {
    grep( paste(cell_type,collapse="|"),colnames(interaction_scores) )
    interaction_scores = interaction_scores[,grep( paste(cell_type,collapse="|"),colnames(interaction_scores) )]
    combined_scores = combined_scores[,grep( paste(cell_type,collapse="|"),colnames(combined_scores) ),]
  }
  
  #Change row names
  rownames(interaction_scores) = gsub("_","-", rownames(interaction_scores))
  colnames(interaction_scores) = gsub("_","-", colnames(interaction_scores))
  
  tests = apply( combined_scores, c(1,2), wilcox.test, alternative = "greater",mu = 0)
  p_values = matrix(do.call(rbind,lapply(tests,function(v){v$p.value})), nrow = dim(tests)[1])
  
  
  FDR_values = matrix(p.adjust( as.vector(p_values), method="BH"), nrow = dim(p_values)[1])
  
  ##### Create score heatmap figure
  # Display only rows/columns with any score greater than threshold
  row_idx = apply(interaction_scores>threshold,1,any)
  col_idx = apply(interaction_scores>threshold,2,any)
  #row_idx = apply( p_values < .02 ,1 , any )
  #col_idx = apply( p_values < .02 ,2 , any )
  interaction_subset = interaction_scores[row_idx,col_idx]
  combined_subset = combined_scores[row_idx,col_idx,]
  p_value_subset = p_values[row_idx,col_idx]
  FDR_subset = FDR_values[row_idx,col_idx]
  
  # Create color map
  cVals = seq( from = 0, to = ceiling(max(interaction_subset)), by = 0.02 )
  color_map = colorRamp2(cVals,col_func(length(cVals)))

  #null_scores = null_scores[row_idx,col_idx,]
  # print figure
  pdf(paste(out_dir,"fig_2",cell_type,"_combined_scores.pdf"),width=11,height=8.5,paper='special') # open pdf 
  hMap = Heatmap(interaction_subset, name ="combined", column_title ="combined",gap = unit(.02,"npc"),
                 row_names_side = "left",
                 na_col = "white", col = color_map, rect_gp = gpar(col = "gray"),
                 cluster_columns = FALSE,cluster_rows = FALSE,
                 show_heatmap_legend = TRUE)
  
  print(hMap) # print heatmap to pdf
  #### Annotate significant interactions
  num_row = dim(interaction_subset)[1]
  num_col = dim(interaction_subset)[2]
  significant_idx = arrayInd( which( FDR_subset < FDR ), dim(interaction_subset)) # indices of significant scores
  sig_row_idx = (num_row + 1) - significant_idx[,1]  # row_idx starts from bottom left
  sig_col_idx = significant_idx[,2]
  # Convert to units for figure
  x_coord = unit((sig_col_idx-0.5)/num_col, "npc")
  y_coord = unit((sig_row_idx-0.5)/num_row, "npc")
  width = unit(1/num_col, "npc")
  height = unit(1/num_row, "npc")
  decorate_heatmap_body("combined",  code = {dec = grid.circle(x_coord, y_coord, r = unit(0.25, "mm"),
                                                          default.units = "npc", name = NULL,
                                                          gp=gpar(col = NULL,alpha=1,lwd=1,lty="solid",fill="black"), draw = TRUE, vp = NULL) })
  
  
  dev.off()
  
  
}