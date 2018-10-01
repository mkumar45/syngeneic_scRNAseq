function trained_classifier = train_classifier( X , training_labels, gene_symbols )
%TRAIN_DECISION_TREE train cell-type classifier
%defined in marker_path
%   INPUTS:
%       X: expression data 
%       training_labels : cell-type labels 
%   OUTPUTS:
%       trained_classifier: classifer with predict function

    % Use only 500 most variable genes and cells with assigned training
    % label
    gene_std = std(X,[],1);
    [~,sort_idx] = sort( gene_std,'descend' );
    gene_subset = @(z) z(:,sort_idx(1:500)); 
    X_training = gene_subset(X);
    
    % Apply PCA to the predictor matrix.
    [pca_coefficients, pca_scores, ~, ~, variance_explained, pca_centers] = pca( X_training, ...
                                                                        'Centered', true,...
                                                                        'Algorithm','svd');

    % Keep enough components to explain the desired amount of variance.
    explained_variance_to_keep = 95/100;
    num_Components = find(cumsum(variance_explained)/sum(variance_explained) >= explained_variance_to_keep, 1);
    pca_coefficients = pca_coefficients(:,1:num_Components);
    predictors = pca_scores(:,1:num_Components);

    % Train classifier
    classification_tree = fitctree(predictors, training_labels, ...
                                  'SplitCriterion', 'gdi', ...
                                   'ClassNames', unique(training_labels), ...
                                   'MaxNumSplits', 100, ...
                                   'Surrogate', 'off');

    % Create the result struct with predict function
    pca_transformation = @(z) bsxfun(@minus, z, pca_centers) * pca_coefficients;
    tree_prediction = @(z) predict(classification_tree, z);
    trained_classifier.predict = @(x) tree_prediction(pca_transformation(gene_subset(x)));    
    
    trained_classifier.predict_confidence = @(z,confidence) predict_confidence(z, trained_classifier.predict, confidence);
    
    % Add additional fields to the result struct
    trained_classifier.PCACenters = pca_centers;
    trained_classifier.PCACoefficients = pca_coefficients;
    trained_classifier.PCAScores = pca_scores;
    trained_classifier.ClassificationTree = classification_tree;
    trained_classifier.VarianceExplained = variance_explained;

    if (exist('gene_symbols','var'))
        trained_classifier.ClassificationGenes = gene_symbols(sort_idx(1:500));
    end
end
