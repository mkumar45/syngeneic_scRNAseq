function [ output_args ] = make_tirosh_phenotype( input_file, cell_type_path, column_path, row_path, output_path )
%MAKE_TIROSH_PHENOTYPE Quantify phenotypes of interest for tirosh dataset
%   Detailed explanation goes here

    % Import data
    X = csvread( input_file );
    cell_type_table = readtable(cell_type_path);
    predicted_cell_type = cell_type_table.predicted_cell_type;
    column_data = readtable( column_path );
    
    %% Percentage of Tregs
    count_tbl = count_cell_types( predicted_cell_type, cellstr( num2str(column_data.tumor_id )) );
    pct_treg_tfac = count_tbl.Treg./ sum( count_tbl{:,{'Treg','Th','CD8'}},2);
    pct_treg_total = count_tbl.Treg ./ sum( count_tbl{:,:},2); 

    %% Percentage of proliferating melanoma cells
    tumor_ids = column_data.tumor_id; 
    tumors = unique(tumor_ids);
    num_tumors = length(tumors);
    mel_idx = strcmp( predicted_cell_type, 'Mel'); % index of melanoma cells

    % Use Ki67 (MKI67) as a marker for proliferation
    row_data = readtable( row_path );
    gene_symbols = row_data.GeneSymbol;
    ki67_idx = ismember( gene_symbols, 'MKI67' );

    % Split melanoma cells in Ki67 +/- populations
    marker_expression = X( mel_idx, ki67_idx );
    k_fold_models = fit_GMM( marker_expression );
    [~, cluster_idx, sort_idx] = select_GMM( marker_expression, k_fold_models, {'MKI67'} );
    proliferating_idx = bsxfun(@eq, cluster_idx, cellfun( @(x) x(1), sort_idx) );

    pct_ki_pos = zeros(num_tumors,1);   
    for t = 1:num_tumors
        tumor_idx =  tumor_ids( mel_idx) == tumors(t) ; % index of patient-specific melanoma cells
        if any ( tumor_idx  )
            pct_ki_pos(t) = sum(tumor_idx & proliferating_idx) ./ sum( tumor_idx );
        else
            pct_ki_pos(t) = NaN;
        end   
    end
    
    %%
    tbl = table( pct_ki_pos, pct_treg_tfac, pct_treg_total ); 
    writetable( tbl, output_path );

    
end

