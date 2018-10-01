function tirosh_analysis(input_file, row_path, column_path, tsne_path, cell_type_path )

%     input_file = 'data/syngeneic_expression.csv';
%     tsne_path = 'results/tirosh_tsne.csv';
%     column_path = 'results/tirosh_processed_column_data.csv';
%     row_path = 'results/tirosh_processed_row_data.csv';
%     RL_path = 'data/receptor_ligand_Ramilowski_modified.csv';
%     cell_type_path = 'results/tirosh_classification/tcell_subsets/cell_type_labels.txt';

    %% Import data
    X = csvread( input_file ); 
    row_data = readtable(row_path);
    column_data = readtable( column_path );
    tsne_coordinates = readtable( tsne_path );

    % Import cell-type predictions and use only cells with assigned labels confidence
    [ cell_type_labels ] = readtable( cell_type_path );
    confident_idx = ~strcmp( cell_type_labels.predicted_cell_type, 'Unassigned' );
    X_subset = X( confident_idx, : ); 
    tsne_subset = tsne_coordinates( confident_idx, : ); 
    
    %% t-SNE plot of T cell subsets
    cMap_ct = color_maps; 
    cMap = cMap_ct([1,12,11,9,13],:);
    tcell_labels = cell_type_labels.predicted_cell_type( confident_idx );
    tcell_labels( ~ismember( tcell_labels, {'CD8','Th','Treg','Tcell'} ) ) = {'Unassigned'};
    
    unique_ct = unique(tcell_labels);
    for ct = 1:length(unique_ct)
        idx = ismember( tcell_labels, unique_ct(ct) );
        scatter( tsne_subset{ idx,1 }, tsne_subset{ idx, 2}, 10,...
                 'MarkerEdgeColor', cMap(ct,:) , 'MarkerFaceColor', cMap(ct,:) );
        
        hold on
    end
    xlabel('t-SNE1'); ylabel('t-SNE2'); axis square; box on;
    legend(unique_ct)
    set(gca,'XTick',[]); set(gca,'YTick',[]);  set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'FontSize',16)
    set(gcf, 'PaperPositionMode', 'auto'); legend('boxoff','Location','NorthEast')    
    saveas(gcf, strcat( 'figures/sup_fig_4_tsne.pdf') )
    
    
    %% Histogram
    cell_type_count = count_cell_types( cell_type_labels.predicted_cell_type,...
                                        cellstr(num2str(column_data.tumor_id)) );
    pct_cell_type = bsxfun( @rdivide, cell_type_count{:,:}, sum(cell_type_count{:,:}, 2));
    figure
    figHand = bar( pct_cell_type, 'stack');
    
    set(gca,'XTickLabels', cell_type_count.Properties.RowNames )
    legend( cell_type_count.Properties.VariableNames );
    box on; axis square; set(gca,'FontSize',16); ylim([0,1])
    xlabel('Patient ID'); ylabel('Cell type percentage');
    set(gca,'YTick',0:0.2:1); set(gca,'YTickLabels',{'0','20','40','60','80','100'})
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gcf, 'PaperPositionMode', 'auto')
    legend('boxoff','Location','NorthEast')
    
    figHand(1).FaceColor = cMap_ct(2,:); 
    figHand(2).FaceColor = cMap_ct(3,:);
    figHand(3).FaceColor = cMap_ct(1,:);
    figHand(4).FaceColor = cMap_ct(6,:);
    figHand(5).FaceColor = cMap_ct(8,:);
    figHand(6).FaceColor = cMap_ct(4,:);
    figHand(7).FaceColor = cMap_ct(10,:);
    figHand(8).FaceColor = cMap_ct(12,:);
    figHand(9).FaceColor = cMap_ct(11,:);
    figHand(10).FaceColor = cMap_ct(9,:);
    figHand(11).FaceColor = cMap_ct(13,:);

    saveas(gcf, strcat( 'figures/sup_fig_4_cell_type_histogram.pdf') )

    %% Predict Treg pecentage
    
    phenotype_tbl = readtable('results/tirosh_phenotype.csv');
    % Import interaction scores
    file_list = dir('results/tirosh_communication/*_interaction_scores.csv');
    num_files = length(file_list);
    
    score_array = cell(num_files,1);
    pair_names = {}; % cell-type pairs

    for n = 1:num_files
        score_array{n} = readtable( strcat('results/tirosh_communication/',file_list(n).name),...
                                    'ReadRowNames',true);
        pair_names = union( pair_names, score_array{n}.Properties.VariableNames, 'sorted' );
    end
    % Add in missing interactions
    score_matrix = zeros( num_files, length(pair_names)*size(score_array{1},1)); 

    for n = 1:num_files
        missing_pairs = setdiff( pair_names, score_array{n}.Properties.VariableNames);
        score_array{n}(:,missing_pairs) = {0};
        score_array{n} = score_array{n}(:,pair_names);
        score_matrix(n,:) = score_array{n}{:,:}(:);
    end
    
    normalized_matrix = bsxfun( @rdivide, bsxfun( @minus, score_matrix, min(score_matrix)), max(score_matrix) - min(score_matrix) );
    
    num_cell_pairs =size(score_array{1},2);
    num_RL = size(score_array{1},1);
    
    interaction_names = strcat( repmat( score_array{1}.Properties.RowNames, num_cell_pairs ,1 )',...
                                '_',repelem( score_array{1}.Properties.VariableNames, num_RL ) );
                            
    idx = sum( normalized_matrix == 0 ) < 10 & ~all(isnan(normalized_matrix)); 
    [B, stats] = lasso( normalized_matrix(:,idx) , phenotype_tbl.pct_treg_tfac,...
                        'CV', 5, 'PredictorNames', interaction_names(idx) );

    plot_lasso(B,stats,'lasso_Treg')
    
    %% Species comparison
    prediction_table = readtable( 'results/tirosh_classification//cell_type_predictions.csv' );
    confidence = 0.95;
    % Use only cells with confident cell-type prediction
    confident_idx = any(prediction_table{:,:} >= confidence,2);
    predicted_cell_type = cell( size( prediction_table, 1 ), 1);
    predicted_cell_type( ~confident_idx) = {'Unassigned'};
    [~, col_idx] = max( prediction_table{confident_idx,:} > confidence,[],2 );    
    predicted_cell_type( confident_idx) = prediction_table.Properties.VariableNames( col_idx );
    cell_types = predicted_cell_type(confident_idx);
    cell_types( strcmp(cell_types, 'Mel')) = {'Tumor'};
    
    gene_symbol = row_data.GeneSymbol; 
    
    RL_list = readtable(RL_path);
    filtered_RL = filter_RL_list( RL_list, gene_symbol );
    
    
    % Score for each tumor
    by_col = 'tumor_id';
    output_path = 'results/tirosh_comparison/';
    column_labels = column_data{ confident_idx , by_col };
    unique_columns = unique( column_labels);

    for col = 1:length(unique_columns)
        column_idx = column_labels == unique_columns(col) ;
        [ interaction_table, null_scores, receptor_table,ligand_table]  = compute_scores( X_subset( column_idx, : ), filtered_RL, cell_types( column_idx) , gene_symbol ); 

        %csvwrite( strcat(output_path,unique_columns{col},'_interaction_scores.csv'), interaction_scores);
        writetable( interaction_table, strcat(output_path,num2str(unique_columns(col)),'_interaction_scores.csv'),...
                    'WriteRowNames',true)
        writetable( receptor_table, strcat(output_path,num2str(unique_columns(col)),'_receptor_expression.csv'),...
                    'WriteRowNames',false)
        writetable( ligand_table, strcat(output_path,num2str(unique_columns(col)),'_ligand_expression.csv'),...
                    'WriteRowNames',false)
        write_null_scores( strcat(output_path,num2str(unique_columns(col)),'_null_scores/'), null_scores );
    end
        


end
