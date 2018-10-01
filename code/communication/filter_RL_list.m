function [ filtered_RL ] = filter_RL_list( RL_list, gene_symbols )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
%% Use only PPIs for which both proteins match GO terms
    % Read in table that contains genes matching specified GO terms
%     gene_list = readtable( file_path, 'ReadRowNames', false,...
%                                       'ReadVariableNames',true,...
%                                       'Delimiter','\t');
%     gene_symbols = unique(gene_list.Symbol,'stable'); % table has multiple entries (per GO term) for each gene
%     [~, go_idx_a] = ismember( interaction_list.gene_symbol_a, gene_symbols ); % can't use intersect -- doesn't return repeitions
%     [~, go_idx_b] = ismember( interaction_list.gene_symbol_b, gene_symbols ); % can't use intersect -- doesn't return repeitions
%     go_idx = go_idx_a & go_idx_b;
%     interaction_list = interaction_list(go_idx,:);
    
    %% Filter out genes with low average expression
%     keep_idx = mean( flexmatrix.data,2 ) > 0.1;
%     sum(keep_idx)
%    select_rows( flexmatrix, keep_idx );

    %% Use only interactions that have measured genes
    receptor_list = RL_list.Receptor_ApprovedSymbol;
    ligand_list = RL_list.Ligand_ApprovedSymbol;
   
    interaction_idx = ismember(receptor_list, gene_symbols ) & ismember(ligand_list, gene_symbols ); % indices of interactions with both genes measured
    filtered_RL = RL_list(interaction_idx,:);
    
    %% Remove duplicate rows and use only literature supported interactions
    interaction_pairs = table(filtered_RL.Receptor_ApprovedSymbol, filtered_RL.Ligand_ApprovedSymbol);
    [~,unique_idx] = unique(interaction_pairs,'rows','stable');
    filtered_RL = filtered_RL( unique_idx,: );
    filtered_RL = filtered_RL( strcmp( filtered_RL.Pair_Evidence, 'literature supported') ,...
                        {'Pair_Name', 'Ligand_ApprovedSymbol', 'Receptor_ApprovedSymbol'});


end

