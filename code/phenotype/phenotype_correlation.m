function [ var, interaction_scores ] = phenotype_correlation( scores_dir, phenotype_path, num_reps, out_dir )
%UNTITLED6 Summary of this function goes here
% INPUTS: 
%   model_flexmat - cell containing score flexmat for each model
    if ~exist('num_reps','var')
        num_reps = 10; % analyze all models by default
    end
      
    file_list = dir( strcat( scores_dir, '*_interaction_scores.csv'));
    receptor_list = dir( strcat( scores_dir, '*_receptor_expression.csv'));
    ligand_list = dir( strcat( scores_dir, '*_ligand_expression.csv'));
    phenotype_table = readtable( phenotype_path );
    num_models = length( file_list );
    
    % Import scores for each model
    interaction_scores = cell( num_models, 1 );
    receptor_expression = cell( num_models, 1 );
    ligand_expression = cell( num_models, 1 );
    
    pair_names = {}; % cell-type pairs
    cell_type_names = {}; % cell-types 
    for n = 1:num_models
        interaction_scores{n} = readtable( strcat( scores_dir, file_list(n).name) , 'ReadRowNames', true );
        pair_names = union( pair_names, interaction_scores{n}.Properties.VariableNames, 'sorted' );
        receptor_expression{n} = readtable( strcat(scores_dir, receptor_list(n).name ) , 'ReadRowNames', false );
        cell_type_names = union( cell_type_names, receptor_expression{n}.Properties.VariableNames, 'sorted' );
        ligand_expression{n} = readtable( strcat(scores_dir,ligand_list(n).name) , 'ReadRowNames', false );
    end
    
    % Add in missing interactoins
    interaction_mat = zeros( num_models, length(pair_names)*size(interaction_scores{1},1)); 
    receptor_mat = zeros( num_models, length(cell_type_names)*size(interaction_scores{1},1)); 
    ligand_mat = zeros( num_models, length(cell_type_names)*size(interaction_scores{1},1) ); 
    
    for n = 1:num_models
        missing_pairs = setdiff( pair_names, interaction_scores{n}.Properties.VariableNames);
        interaction_scores{n}(:,missing_pairs) = {0};
        interaction_scores{n} = interaction_scores{n}(:,pair_names);
        interaction_mat(n,:) = interaction_scores{n}{:,:}(:);
        
        missing_ct = setdiff( cell_type_names, receptor_expression{n}.Properties.VariableNames);
        receptor_expression{n}(:,missing_ct) = {0};
        receptor_expression{n} = receptor_expression{n}(:,cell_type_names);
        receptor_mat(n,:) = receptor_expression{n}{:,:}(:);
        ligand_expression{n}(:,missing_ct) = {0};
        ligand_expression{n} = ligand_expression{n}(:,cell_type_names);
        ligand_mat(n,:) = ligand_expression{n}{:,:}(:);  
    end
    
    [scores_correlation, scores_p_value] = corr( interaction_mat, phenotype_table{:,:} , 'type', 'spearman', 'rows', 'pairwise' );
    [receptor_correlation, receptor_p_value] = corr( receptor_mat, phenotype_table{:,:} , 'type', 'spearman', 'rows', 'pairwise' );
    [ligand_correlation, ligand_p_value] = corr( ligand_mat, phenotype_table{:,:} , 'type', 'spearman', 'rows', 'pairwise' );
    
    % Write files
    interaction_names = interaction_scores{1}.Properties.RowNames;
    for p = 1:size(phenotype_table,2)
        phenotype_name = phenotype_table.Properties.VariableNames{p};
        % Interaction scores
        tbl = array2table( reshape( scores_correlation(:,p), [], length(pair_names)),...
                           'VariableNames', pair_names, ...
                           'RowNames', interaction_names );
        writetable( tbl, strcat( out_dir, 'score_correlation_',phenotype_name,'.csv'), 'WriteRowNames',true )
        tbl = array2table( reshape( scores_p_value(:,p), [], length(pair_names)),...
                           'VariableNames', pair_names, ...
                           'RowNames', interaction_names );
        writetable( tbl, strcat( out_dir, 'score_p_value_',phenotype_name,'.csv'), 'WriteRowNames',true )
        % Receptor Expression
        tbl = array2table( reshape( receptor_correlation(:,p), [], length(cell_type_names)),...
                           'VariableNames', cell_type_names, ...
                           'RowNames', interaction_names );
        writetable( tbl, strcat( out_dir, 'receptor_correlation_',phenotype_name,'.csv'), 'WriteRowNames',true )
        tbl = array2table( reshape( receptor_p_value(:,p), [], length(cell_type_names)),...
                           'VariableNames', cell_type_names, ...
                           'RowNames', interaction_names );
        writetable( tbl, strcat( out_dir, 'receptor_p_value_',phenotype_name,'.csv'), 'WriteRowNames',true )
        % Ligand Expression
        tbl = array2table( reshape( ligand_correlation(:,p), [], length(cell_type_names)),...
             'VariableNames', cell_type_names, ...
             'RowNames', interaction_names );
        writetable( tbl, strcat( out_dir, 'ligand_correlation_',phenotype_name,'.csv'), 'WriteRowNames',true )
        tbl = array2table( reshape( ligand_p_value(:,p), [], length(cell_type_names)),...
                           'VariableNames', cell_type_names, ...
                           'RowNames', interaction_names );
        writetable( tbl, strcat( out_dir, 'ligand_p_value_',phenotype_name,'.csv'), 'WriteRowNames',true )
    end
    
end
    
