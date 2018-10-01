function [ null_corr ] = compute_null_scores( scores, phenotype_mat, num_reps )
%UNTITLED6 Summary of this function goes here
% INPUTS: 
%   model_flexmat - cell containing score flexmat for each model

    num_rows = size(phenotype_mat,1);
    
    num_perms = factorial( num_rows ); % number of possible phenotype permutations
    num_interactions = size( scores , 2); 
    num_phenotypes = size(phenotype_mat,2);
    
    null_corr = zeros( num_interactions, num_phenotypes , num_reps );

    
    perm_idx = perms( 1: num_rows ); % inidces of all possible permutations
    rand_idx = perm_idx(randperm( num_perms, num_reps ),: ); % pick random subsest of permutatoins
    
    for p = 1:num_phenotypes
        phenotype_vec = phenotype_mat(:,p);
        size(phenotype_vec(rand_idx'))
        size( scores) 
        [rho, ~] = corr( scores, phenotype_vec(rand_idx') , 'type', 'spearman', 'rows', 'pairwise' ); 
        null_corr(:,p,:) = rho; 
    end


end