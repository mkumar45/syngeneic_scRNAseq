function [ cell_type_accuracy, cell_type_count ] = training_cross_validatation( X ,labels,k_fold )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Perform cross-validation

    unique_cell_types = unique(labels);
    num_cell_types = length(unique_cell_types);
     
    if(~exist('k_fold','var'))
        k_fold = 5;
    end
    
    cvp = cvpartition(labels, 'KFold', k_fold);
    
    cell_type_accuracy = zeros(k_fold, num_cell_types); 
    cell_type_count = zeros(k_fold, num_cell_types);
    fold_accuracy = zeros(k_fold,1);

    for k = 1:k_fold
        % Split into training/testing data set
        X_training = X(cvp.training(k),:);
        X_test = X(cvp.test(k),:);
        training_labels = labels(cvp.training(k));
        test_labels = labels(cvp.test(k));

        cell_type_count(k, :) =  cellfun( @(x) sum(ismember( training_labels ,x)), unique_cell_types); 
                
        % Train classifier
        fold_classifier = train_classifier( X_training, training_labels );
        
        % Compute validation predictions and scores
        [fold_predictions, confidence_idx] = fold_classifier.predict_confidence( X_test, 0.95 ); 

        correct_idx = strcmp(fold_predictions(confidence_idx), test_labels(confidence_idx));
        cell_type_accuracy(k,:) = compute_cell_type_accuracy( unique_cell_types,...
                                                              test_labels(confidence_idx),...
                                                              fold_predictions(confidence_idx) );

                                                          % Assess accuracy of high-confidence predictions
        num_confident = sum(confidence_idx); 
        test_size = length(test_labels);
        fold_accuracy(k) = sum(correct_idx)./ num_confident ;
        
        fprintf('%f of %i high confidence predictions correct. ',fold_accuracy(k),num_confident);
        fprintf('%i out of %i low confidence predictions\n', test_size - num_confident,test_size);
    end
    
end

function pct_correct = compute_cell_type_accuracy( unique_cell_types, labels, predictions )
    pct_correct= NaN(1,length(unique_cell_types));
    for ct = 1:length(unique_cell_types)
        
        ct_idx = strcmp( labels, unique_cell_types{ct} );
        correct = strcmp( predictions( ct_idx), unique_cell_types{ct} );
        pct_correct(ct) = sum(correct)./length(correct);
    end
end

