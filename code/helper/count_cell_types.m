function [ count_tbl ] = count_cell_types( predicted_cell_type, col_labels )
%[ count_tbl ] = count_cell_types( predicted_cell_type, col_labels )
%   Count number of each cell type based on unique labels in col_labels
%   INPUTS: predicted_cell_type = array of predicted cell-types
%           col_labels: corresponding labels to subset for counting
%   OUTPUTS: count_tbl - table of cell-type counts for each label in
%                        col_labels
    
    unique_labels = unique( col_labels );
    cell_types = unique(predicted_cell_type);
    count_mat = zeros( length(unique_labels), length(cell_types) );

    for n = 1:length(unique_labels)
        label = unique_labels(n);
        for ct = 1:length(cell_types)
            cell_type = cell_types(ct);
            count_mat(n,ct) = sum( strcmp(col_labels,label) & strcmp( predicted_cell_type, cell_type ) );
        end
    end

    count_tbl = array2table(count_mat);
    count_tbl.Properties.VariableNames = cell_types;
    count_tbl.Properties.RowNames = unique_labels;


end
