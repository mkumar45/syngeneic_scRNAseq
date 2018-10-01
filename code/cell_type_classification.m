function cell_type_classification( input_file, row_path, marker_path, output_path )
%CELL_TYPE_CLASSIFICATION Classify cell-types defined in marker_table using
%gene expression in input_file
% input_file: .csv file containing containing scRNA-seq expression data.
%              Rows are observations (single-cells) and columns are variables (genes)
% marker_table: .csv file defining cell-types of interest by marker genes
% output_file: path to .csv file to write output
%   Detailed explanation goes here
    
    %% Read in data
    X = csvread( input_file );
    row_data = readtable( row_path );
    gene_symbols = row_data.GeneSymbol;
    
    %% 
    [ X_training, training_labels, mixture_models,idx, parameter_tbl ] = create_training_data( X, gene_symbols, marker_path );
    [ classifier ] = train_classifier( X_training, training_labels, gene_symbols );
    [ cell_type_accuracy, cell_type_count ] = training_cross_validatation( X_training ,training_labels,5 );
    [~, ~, prediction_scores ] = classifier.predict_confidence( X, 0.95 );
    
    
    %% Write results to file
    % Write predicted cell-types and scores to file
    writetable( array2table(prediction_scores,'VariableNames',unique(training_labels)) , strcat(output_path,'cell_type_predictions.csv') );
    % Save gaussian mixture models to .mat file
    save( strcat(output_path,'mixture_models.mat'), 'mixture_models')
    % Write classifier information
    save( strcat(output_path,'classifier.mat'), 'classifier')

    % Write cross-validation metrics
    save( strcat( output_path, 'CV_metrics.mat' ), 'cell_type_accuracy',...
                                                   'cell_type_count' )
                                               
    % GMM paramters
    writetable( parameter_tbl, 'results/syngeneic_classification/GMM_parameters.csv', 'WriteRowNames', true );

    % Confusion matrix
    [ predicted_labels, idx ] = classifier.predict_confidence( X_training, 0.95 )
    cell_types = unique( predicted_labels(idx) ); 
    confusion_mat = array2table( confusionmat( training_labels( idx), predicted_labels(idx), 'ORDER', cell_types ) ,...
                                 'VariableNames', cell_types, 'RowNames', cell_types );
    writetable( confusion_mat , strcat(output_path,'confusion_matrix.csv'),...
                'WriteRowNames',true,'WriteVariableNames',true);

end