function tirosh_communication( input_file, row_path, column_path, cell_type_path, RL_path, by_col, output_path )
%TIROSH_COMMUNICATION Summary of this function goes here
%   Detailed explanation goes here
%     input_file = 'data/syngeneic_expression.csv';
%     row_path = 'data/syngeneic_row_data.csv';
%     column_path = 'data/syngeneic_column_data.csv';
%     cell_type_path = 'results/syngeneic_classification/cell_type_predictions.csv';
%     RL_path = 'data/mouse_receptor_ligand.csv';
%     by_col = 'Model';

    X = csvread( input_file );  
    row_data = readtable( row_path );
    column_data = readtable( column_path );
    gene_symbols = row_data.GeneSymbol;

    % Import cell-type predictions and use only cells with assigned labels confidence
    [ cell_type_labels ] = readtable( cell_type_path );
    confident_idx = ~strcmp( cell_type_labels.predicted_cell_type, 'Unassigned' );
    X_subset = X( confident_idx, : ); 
    cell_type_labels = cell_type_labels.predicted_cell_type( confident_idx ); 
    
    RL_list = readtable(RL_path);
    filtered_RL = filter_RL_list( RL_list, gene_symbols );
    if ~exist( 'by_col', 'var' )
        compute_scores( X_subset, filtered_RL, cell_type_labels, gene_symbols )
    else
        if sum( strcmp( by_col, column_data.Properties.VariableNames ) ) == 1
            column_labels = column_data{ confident_idx , by_col };
            unique_columns = unique( column_labels);
            
            for col = 1:length(unique_columns)
                column_idx = column_labels == unique_columns(col) ;
                [ interaction_table, null_scores, receptor_table,ligand_table]  = compute_scores( X_subset( column_idx, : ), filtered_RL, cell_type_labels(column_idx), gene_symbols ); 
                
                %csvwrite( strcat(output_path,unique_columns{col},'_interaction_scores.csv'), interaction_scores);
                writetable( interaction_table, strcat(output_path,num2str(unique_columns(col)),'_interaction_scores.csv'),...
                            'WriteRowNames',true)
                writetable( receptor_table, strcat(output_path,num2str(unique_columns(col)),'_receptor_expression.csv'),...
                            'WriteRowNames',false)
                writetable( ligand_table, strcat(output_path,num2str(unique_columns(col)),'_ligand_expression.csv'),...
                            'WriteRowNames',false)
                write_null_scores( strcat(output_path,num2str(unique_columns(col)),'_null_scores/'), null_scores );
            end
        else
            error('Input "by_col" does not match any variable in column_data or is not unique')
        end            

    end
end




