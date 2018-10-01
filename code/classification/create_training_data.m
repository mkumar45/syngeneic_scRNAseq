function [ X_training, training_cell_types, mixture_models, training_idx, parameter_table ] = create_training_data( X, gene_symbols, marker_path )
%CREATE_TRAINING_DATA create training data using cell-types and markers
%defined in marker_path
%   INPUTS:
%       X: expression data 
%       gene_symbols: gene_symbols for each row in X
%       marker_path: filepath to user definied marker path
%   OUTPUTS:
%       X_training: data to use for training classifier
%       cell_types: cell array containing cell-type training labels. 

      
    num_cells = size( X,1 );

    % Import marker table
    marker_table = read_marker_table( marker_path );
    marker_genes = marker_table.Properties.VariableNames;  
    cell_type_names = marker_table.Properties.RowNames;
    num_cell_types = length(cell_type_names);
     
    % Log-normalize
    X_normalized = log2( 1 + X );

    % Find indices of cell-type markers
    [~, marker_table_idx, flexmat_marker_idx] = intersect( marker_genes, gene_symbols, 'stable');
    marker_expression = X_normalized(:,flexmat_marker_idx);
    marker_table = marker_table(:,marker_table_idx); % Use only measured markers

    k_fold_models = fit_GMM( marker_expression );
    [mixture_models, cluster_idx, sort_idx, parameter_table] = select_GMM( marker_expression, k_fold_models, marker_genes );
    
    % Find cells that match the specified marker profiles
    is_cell_type = zeros(num_cell_types,num_cells);    
    for ct = 1:length(cell_type_names)
        cell_type_profile = marker_table{ct,:}; % marker expression pattern for cell-type
        idx = ~cellfun(@isempty,cell_type_profile);  % index of markers for cell-type 

        if any(idx)
            and_idx = strcmp(cell_type_profile,'1');
            not_idx = strcmp(cell_type_profile,'0');
            
            % Find cells matches marker expression profile              
            all_markers = all( bsxfun(@eq, cluster_idx(:,and_idx), cellfun( @(x) x(1), sort_idx(and_idx)) ), 2 );
            not_markers = all( bsxfun(@eq, cluster_idx(:,not_idx), cellfun( @(x) x(end), sort_idx(not_idx)) ), 2 );
            
            
            if any(and_idx) && any(not_idx)
                is_cell_type(ct,:) = all_markers & not_markers;
            elseif any(and_idx)
                    is_cell_type(ct,:) = all_markers;
            elseif any(not_idx)
                is_cell_type(ct,:) = not_markers;
            else
               error('No markers') 
            end
                 
        end
    end
    
    % Select identified cells         
    training_idx = sum( is_cell_type, 1) == 1;  % only use cells that only have one cell-type above threshold

    % Find cell-type
    [~,max_idx] = max(is_cell_type,[],1);
    cell_types = [cell_type_names(max_idx)];
    
    training_cell_types = cell_types(training_idx);
    X_training = X( training_idx, : );
    
end