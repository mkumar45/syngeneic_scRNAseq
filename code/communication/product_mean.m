function [ interaction_scores ] = product_mean( receptor_expression, ligand_expression )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    num_cell_types = size(receptor_expression,2);
    num_interactions = size(receptor_expression,1);
    num_pairs = num_cell_types.^2 ; % all cell-type pairs + self-interactions
    interaction_scores = zeros(num_interactions,num_pairs);


    receptor_mean_expr = cellfun(@mean,receptor_expression);
    ligand_mean_expr = cellfun(@mean,ligand_expression);

    for n = 1:num_interactions
        score_matrix = receptor_mean_expr(n,:)'*ligand_mean_expr(n,:);
        interaction_scores(n,:) = reshape(score_matrix,1,num_pairs);
    end


end

