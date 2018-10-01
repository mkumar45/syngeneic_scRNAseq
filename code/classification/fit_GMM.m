function [ GMMs ] = fit_GMM( X, num_components, k_fold )
%FIT_GMM Fit GMMS for each row in X. Fits model for each value of num_components. 
% Splits data in X into k_fold partitions and fits k-models for each gene
% and number of componenets. 
%   INPUTS:
%       X: expression data (used for clustering after GMM is selected
%       num_components: vector containing number of compoenents to use for each GMM
%       k_fold: number of partitions to split data for fitting
%   OUTPUTS:
%       GMMs: fit gaussian mixture models
    
    if(~exist('num_classes','var') || isempty(num_components))
        num_components = 5;
    end
    if(~exist('k_fold','var') || isempty(k_fold))
        k_fold = 5;
    end

    num_cells = size(X,1);
    num_genes = size(X,2);

    GMMs = cell(num_genes,num_components);   
    idx = crossvalind('Kfold',num_cells,k_fold);
    options = statset('MaxIter',500);

    for g = 1:num_genes
        for n = 1:num_components
            for k = 1:k_fold
                GMMs{g,n,k} = fitgmdist(X( idx==k ,g ),n,'Options',options,'CovarianceType','diagonal','RegularizationValue',0.01);
            end

        end
    end
end

