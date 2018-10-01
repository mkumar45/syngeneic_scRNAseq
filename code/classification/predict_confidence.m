function [ prediction_labels, high_confidence_idx, prediction_scores ] = predict_confidence( X, prediction_function, confidence )
%PREDICT_CONFIDENCE Predict cell-type labels and return only
%high-confidence prediction greater than specified threshold
%   INPUTS:
%       X - input data for classifier
%       prediction_function - prediction function returned by training
%       classifier
%       confidence - confidence threshold for predictions

    % Compute validation predictions and scores
    [ prediction_labels , prediction_scores ] = prediction_function( X );
    high_confidence_idx = any(prediction_scores >= confidence,2);  
    prediction_labels(~high_confidence_idx) = {'Unassigned'};

end

