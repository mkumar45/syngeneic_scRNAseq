function [ output_args ] = tcell_classification( input_file, cell_type_path, row_path, marker_path, output_path, confidence )
%TCELL_CLASSIFICATION Summary of this function goes here
%   Detailed explanation goes here

    X = csvread( input_file );
    prediction_table = readtable( cell_type_path );
    
    % Use only cells with confident cell-type prediction
    confident_idx = any(prediction_table{:,:} >= confidence,2);
    predicted_cell_type = cell( size( prediction_table, 1 ), 1);
    predicted_cell_type( ~confident_idx) = {'Unassigned'};
    [~, col_idx] = max( prediction_table{confident_idx,:} > confidence,[],2 );    
    predicted_cell_type( confident_idx) = prediction_table.Properties.VariableNames( col_idx );
    
    
    % Pass only cells confidently classified as T-cells
    tcell_idx = strcmp(predicted_cell_type,'Tcell');
    tcell_expression_path = strcat(output_path,'tcell_expression.csv');
    csvwrite( tcell_expression_path, X( tcell_idx, : ) ); 
      
    % Train classifier to predict T-cell subsets
    cell_type_classification( tcell_expression_path, row_path, marker_path, output_path )
    
    % Import T-cell predictions
    tcell_predictions = readtable( strcat(output_path, 'cell_type_predictions.csv' ) );
    
    % Use only cells with confident cell-type prediction
    confident_idx = any(tcell_predictions{:,:} >= confidence,2);
    predicted_tcell_subtype = cell( size( tcell_predictions, 1 ), 1);
    predicted_tcell_subtype( ~confident_idx) = {'Tcell'};
    [~, col_idx] = max( tcell_predictions{confident_idx,:} > confidence,[],2 );    
    predicted_tcell_subtype( confident_idx) = tcell_predictions.Properties.VariableNames( col_idx );
    
    %% 
    predicted_cell_type( tcell_idx) = predicted_tcell_subtype; 
    writetable( table(predicted_cell_type), strcat( output_path, 'cell_type_labels' ) ); 
    
end

