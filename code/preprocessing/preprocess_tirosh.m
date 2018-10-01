function preprocess_tirosh( input_file,  threshold )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % Import data
    tirosh_tbl = readtable(input_file,...
                           'ReadVariableNames',true);    
                        
    X = transpose( tirosh_tbl{ 4:end, 2:end } );  % single cells as rows, genes as columns
    row_data = table( tirosh_tbl{ 4:end,1 }, 'VariableNames', {'GeneSymbol'}'); 
    gene_symbols = row_data.GeneSymbol;
    
    column_data = array2table( transpose(tirosh_tbl{ 1:3, 2:end }) ,...
                              'VariableNames', {'tumor_id','malignant_flag','non_malignant_id'}); 
    
    hk_table = readtable('data/HK_genes.txt','Delimiter','\t'); 

    % Remove cells with low number of genes detected
    idx = sum( X > 0 ,2 ) > threshold;  
    X_subset = X(idx,:);
    column_subset = column_data( idx, : ); 
    
    % Normalize to houes-keeping expression
    human_hk_symbol = cellfun(@strtrim, hk_table{:,1} ,'UniformOutput',false); % Remove trailing white space from gene symbols  
    hk_idx = ismember( gene_symbols, human_hk_symbol ); 
    hk_expression = mean( X_subset( :, hk_idx ), 2 );
    
    % Divide by house keeping expression and log-normalize
    X_normalized = log2( 1 + bsxfun( @rdivide, X_subset, hk_expression ) );
    
    % Write files
    csvwrite( 'results/tirosh_processed_expression.csv', X_normalized );
    writetable( column_subset, 'results/tirosh_processed_column_data.csv' );
    writetable( row_data, 'results/tirosh_processed_row_data.csv' );


end

