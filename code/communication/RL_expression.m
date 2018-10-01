function [ receptor_expression, ligand_expression ] = RL_expression( X, RL_list, predicted_cell_types, gene_symbols )
% [ receptor_expression, ligand_expression ] = RL_expression( X, RL_list, predicted_cell_type, gene_symbols )
%   Returns cell-type specific expression for each gene listed in a list of
%   protein-protein interactions. 
    
        unique_cell_types = unique( predicted_cell_types ); 
        num_cell_types = length(unique_cell_types);
        num_interactions = size(RL_list,1);

        %% Indices of genes in flexmatrix
        receptor_idx = find_indices( RL_list.Receptor_ApprovedSymbol, gene_symbols);
        ligand_idx = find_indices( RL_list.Ligand_ApprovedSymbol, gene_symbols);

        %% Define variables
        receptor_expression = cell(num_interactions, num_cell_types ); 
        ligand_expression = cell(num_interactions, num_cell_types ); 

        %% Calculate cell-type specific expression
        for ct = 1:num_cell_types
            cell_type_idx = strcmp( predicted_cell_types, unique_cell_types{ct} );
            
            if ~isempty( cell_type_idx )
                receptor_expression(:,ct) = num2cell(X( cell_type_idx , receptor_idx ),1);
                ligand_expression(:,ct) = num2cell(X( cell_type_idx , ligand_idx ),1);
            end
        end
    
end


% function [ gene_a_expression, gene_b_expression ] = ppi_expression_cluster(ppi_list, expression_data, predicted_cell_type, gene_symbols, ct_to_return)
%         
%     if isempty( ct_to_return ) || ~exist('ct_to_return','var')
%         cell_types = unique(predicted_cell_type);
%     else
%         cell_types = ct_to_return;
%     end
%         
%     num_cell_types = length(cell_types);
%     num_interactions = size(ppi_list,1);
% 
%     %% Indices of genes in flexmatrix
%     flex_idx_a = find_indices( ppi_list.gene_symbol_a, gene_symbols);
%     flex_idx_b = find_indices( ppi_list.gene_symbol_b, gene_symbols);
% 
%     %% Define variables
%     gene_a_expression = cell(num_interactions, num_cell_types ); 
%     gene_b_expression = cell(num_interactions, num_cell_types ); 
% 
%     %% Calculate cell-type specific expression
%     for ct = 1:num_cell_types
%         
%         cell_type_idx = strcmp( predicted_cell_type, cell_types{ct} );
%         cell_type_data = expression_data( :, cell_type_idx );
%         if ~isempty( cell_type_data )
%             gene_a_expression(:,ct) = num2cell(full(cell_type_data( flex_idx_a, : )),2);
%             gene_b_expression(:,ct) = num2cell(full(cell_type_data( flex_idx_b, : )),2);
%         end
%     end
% end
