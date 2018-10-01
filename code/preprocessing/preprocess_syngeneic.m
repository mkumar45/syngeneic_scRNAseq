function preprocess_syngeneic( input_file, row_path, column_path, threshold, output_file )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % Import data
    X = csvread( input_file );
    row_data = readtable( row_path ); 
    column_data = readtable( column_path ); 
    gene_symbols = row_data.GeneSymbol;
    hk_table = readtable('data/HK_genes.txt','Delimiter','\t'); 

    % Remove cells with low number of genes detected
    idx = sum( X > 0 ,2 ) > threshold;  
    X_subset = X(idx,:);
    column_subset = column_data( idx, : ); 
    % Normalize to houes-keeping expression
    human_hk_symbol = cellfun(@strtrim, hk_table{:,1} ,'UniformOutput',false); % Remove trailing white space from gene symbols
    mouse_hk_symbol = find_homologs( human_hk_symbol, 'human' );
    unique_func = @(x) length(x) == 1;   
    unique_mouse_hk_idx = cellfun( unique_func, mouse_hk_symbol);
    unique_mouse_hk = mouse_hk_symbol( unique_mouse_hk_idx );
    unravel_func = @(x) x{1};
    unique_mouse_hk = cellfun( unravel_func, unique_mouse_hk, 'UniformOutput',false );
    
    hk_idx = ismember( gene_symbols, unique_mouse_hk ); 
    hk_expression = mean( X_subset( :, hk_idx ), 2 );
    % Divide by house keeping expression and log-normalize
    X_normalized = log2( 1 + bsxfun( @rdivide, X_subset, hk_expression ) );
    csvwrite( output_file, X_normalized );
    writetable( column_subset, 'results/syngeneic_processed_column_data.csv' );


end

