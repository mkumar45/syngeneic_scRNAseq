function [ cell_type_labels, return_idx, predicted_cell_type ] = read_cell_type_predictions( cell_type_path, column_path, confidence )
% [ cell_type_labels, return_idx,predicted_cell_type ] = read_cell_type_predictions( cell_type_path, column_path, confidence ) 
%   Import cell type classification results
%   INPUTS:
%       - cell_type_path: path to output results of cell-type
%       classification
%       - column_path: path to column data
%       - confidence: only assign cell-type labels to cells with
%       prediction confidence

    prediction_table = readtable( cell_type_path );
    column_data = readtable( column_path );     
    
    % Use only cells with confident cell-type prediction
    confident_idx = any(prediction_table{:,:} >= confidence,2);
    predicted_cell_type = cell( size( prediction_table, 1 ), 1);
    predicted_cell_type( ~confident_idx) = {'Unassigned'};
    [~, col_idx] = max( prediction_table{confident_idx,:} > confidence,[],2 );    
    predicted_cell_type( confident_idx) = prediction_table.Properties.VariableNames( col_idx );
    
    % Remove predicted tumor cells that don't match model
    mouse_model = column_data{:,'Model'};
    tumor_idx = ismember(predicted_cell_type,unique(mouse_model));
    correct_model = strcmp( predicted_cell_type, mouse_model);
    
    return_idx = (~tumor_idx | (tumor_idx & correct_model)) & confident_idx;

    % Change tumor-model names to "tumor" for model comparison
    cell_type_labels = predicted_cell_type;
    cell_type_labels( tumor_idx ) = {'Tumor'};

end
