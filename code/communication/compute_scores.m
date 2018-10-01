function [ interaction_table, null_scores, receptor_table, ligand_table ] = compute_scores( X, RL_list, cell_type_labels, gene_symbols, num_null_reps )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    if ~exist('num_null_reps','var') || isempty(num_reps)
        num_null_reps = 100;
    end     

    unique_cell_types = unique(cell_type_labels); 
    num_cell_types = length( unique_cell_types );
    num_interactions = size(RL_list,1);
    num_pairs = num_cell_types.^2 ; % all cell-type pairs + self-interactions
    num_cells = size(X,1);
    
    [ receptor_expression, ligand_expression ] = RL_expression( X , RL_list, cell_type_labels, gene_symbols );
    interaction_scores = product_mean( receptor_expression, ligand_expression );

    %% Calculate null scores
    null_scores = zeros(num_interactions,num_pairs, num_null_reps);
    parfor n = 1:num_null_reps
        rand_idx = randperm( num_cells );
        [ receptor_null_expression, ligand_null_expression ] = RL_expression( X , RL_list, cell_type_labels(rand_idx), gene_symbols );
        null_scores(:,:,n) = product_mean( receptor_null_expression, ligand_null_expression );
    end     
    
    cell_type_pairs = strcat( unique_cell_types(repelem(1:num_cell_types,num_cell_types)),...
                               '_', unique_cell_types( repmat(1:num_cell_types,1,num_cell_types)) );
   interaction_table = array2table( interaction_scores, 'RowNames', strcat(RL_list.Ligand_ApprovedSymbol,'_', RL_list.Receptor_ApprovedSymbol),...
                                    'VariableNames', cell_type_pairs' );
   
   receptor_table = array2table( cellfun(@mean, receptor_expression), 'VariableNames', unique_cell_types' );
   ligand_table = array2table( cellfun(@mean, ligand_expression), 'VariableNames', unique_cell_types' );
   
   
end