function [ mappedX ] = run_tsne( input_file, output_file )
%RUN_TSNE Wrapper function that runs fast_tsne.m (implementation of
%bh_tsne)
% input_file: .csv file containing containing scRNA-seq expression data.
%              Rows are observations (single-cells) and columns are variables (genes)
% output_file: path to .csv file to write output
%   Detailed explanation goes here
    
    X = csvread( input_file );
    fprintf('Loaded Data. %i Rows by %i columns\n',size(X,1),size(X,2));

    gene_std = std(X,[],1);
    [~,sort_idx] = sort( gene_std,'descend' );
    X_subset = X(:,sort_idx(1:500)); 
    
    fprintf('Running t-SNE\n')
    mappedX = fast_tsne(X_subset);
    fprintf('Done t-SNE\n')
    csvwrite(output_file,mappedX)
    fprintf('Output file written\n')
end

