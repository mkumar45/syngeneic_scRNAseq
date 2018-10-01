function [ tcell_table ] = compute_t_cell_pct( cell_type_path, column_path )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

    column_data = readtable( column_path );
    [ cell_type_labels, confident_idx ] = read_cell_type_predictions( cell_type_path, column_path, 0.95 );

    protocol = strrep( column_data.Protocol, 'overnight', 'fresh' ); % overnight samples are unsorted
    unsort_idx = strcmp( protocol, 'fresh' );
    
   
   unique_models = unique(column_data.Model);
   num_models = length(unique_models);
   pct_tcell = zeros(num_models,1);

   
   
   for m = 1:num_models
       model_idx = strcmp( column_data.Model( unsort_idx & confident_idx) , unique_models{m} ); 
       pct_tcell(m) = sum( strcmp( cell_type_labels( model_idx ), 'Tcell') )./ sum(model_idx);
   end
   tcell_table = table( pct_tcell, 'RowNames', unique_models);
       
    
end

