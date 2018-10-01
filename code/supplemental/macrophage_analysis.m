function macrophage_analysis(input_file, row_path, column_path, tsne_path, cell_type_path )
% input_file = 'results/syngeneic_processed_expression.csv'; row_path = 'data/syngeneic_row_data.csv'; column_path = 'results/syngeneic_processed_column_data.csv'; cell_type_path = 'results/syngeneic_classification/cell_type_predictions.csv';
    
    % Import data
    X = csvread( input_file ); 
    row_data = readtable( row_path );
    column_data = readtable( column_path );
    tsne_coordiantes = readtable('results/syngeneic_tsne.csv');
    [ cMap_ct, cMap_model, ~ ] = color_maps();
% Import predicted cell-types
    [ cell_type_labels, confident_idx ] = read_cell_type_predictions( cell_type_path, column_path, 0.95 );
    
    %% DBSCAN
    cluster_idx = dbscan( tsne_coordiantes{:,:}, 2.5, 25); 
    cluster_label = cellstr( strcat( 'Cluster',{' '},num2str( cluster_idx ) ) );
    cluster_label( strcmp( cluster_label, 'Cluster  0' ) ) = {'Unassigned'};
    
    gscatter( tsne_coordiantes{:,1}, tsne_coordiantes{:,2}, cluster_label, ...
              cMap_ct([1,2,3,4,5,6,13,8,9,10,11,12,7,14],:),'.',15 )
    legend( unique(cluster_label,'stable') )
    xlabel('t-SNE1'); ylabel('t-SNE2'); axis square; box on;
    set(gca,'XTick',[]); set(gca,'YTick',[]);  set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'FontSize',16)
    set(gcf, 'PaperPositionMode', 'auto'); legend('boxoff','Location','NorthEast')  
    saveas(gcf, strcat( 'figures/', 'sup_fig_2_tumor_dbscan.pdf') )    
    %% Macrophge analysis
    macrophage_idx = strcmp(cell_type_labels, 'M');
    
    X_macrophage = X( macrophage_idx, : );
    column_data_macrophage = column_data( macrophage_idx,:);
    tsne_macrophage = tsne_coordiantes{macrophage_idx,:};
    
    mouse_strain = column_data_macrophage.Model; 
    mouse_strain( strcmp(mouse_strain, 'B16F10') ) = {'C57BL/6'};
    mouse_strain( strcmp(mouse_strain, 'LL2') ) = {'C57BL/6'};
    mouse_strain( strcmp(mouse_strain, 'MC38') ) = {'C57BL/6'};
    mouse_strain( strcmp(mouse_strain, 'CT26') ) = {'BALB/c'};
    mouse_strain( strcmp(mouse_strain, 'EMT6') ) = {'BALB/c'};
    mouse_strain( strcmp(mouse_strain, 'SA1N') ) = {'aj'};
    
    gscatter(tsne_macrophage(:,1), tsne_macrophage(:,2), mouse_strain, cMap_ct([1,4,11],:),'.',15  )
    xlabel('t-SNE1'); ylabel('t-SNE2'); axis square; box on;
    set(gca,'XTick',[]); set(gca,'YTick',[]);  set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'FontSize',16)
    set(gcf, 'PaperPositionMode', 'auto'); legend('boxoff','Location','NorthEast')  
    saveas(gcf, strcat( 'figures/', 'sup_fig_1_macrophage_model.pdf') )    

    %% DBSCAN
    cluster_idx = dbscan( tsne_macrophage, 2, 30); 

    %% Plots
    cMap_ct(14,:) = [ 255,204,10]./256;
    cluster_label = cellstr( strcat( 'Cluster',{' '},num2str( cluster_idx ) ) )
    cluster_label( strcmp( cluster_label, 'Cluster 0' ) ) = {'Unassigned'};
    figure
    subplot(121)
    gscatter(tsne_macrophage(:,1), tsne_macrophage(:,2), cluster_label, cMap_ct( [2,13,3,6,8,10,12,14],: ),'.',15 )
    xlabel('t-SNE1'); ylabel('t-SNE2'); axis square; box on;
    set(gca,'XTick',[]); set(gca,'YTick',[]);  set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'FontSize',16)
    set(gcf, 'PaperPositionMode', 'auto'); legend('boxoff','Location','NorthEast')  
    legend(unique(cluster_label,'stable'))
    %saveas(gcf, strcat( out_path, 'fig_1_tsne_model.pdf') )    subplot(122)
    subplot(122)
    gscatter(tsne_macrophage(:,1), tsne_macrophage(:,2), column_data_macrophage.Model,cMap_model,'.',15 )
    xlabel('t-SNE1'); ylabel('t-SNE2'); axis square; box on;
    set(gca,'XTick',[]); set(gca,'YTick',[]);  set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    set(gca,'FontSize',16)
    set(gcf, 'PaperPositionMode', 'auto'); legend('boxoff','Location','NorthEast')  
    legend(unique(column_data_macrophage.Model,'stable'))
    saveas(gcf, strcat( 'figures/', 'sup_fig_1_macrophage.pdf') )    
       
    %%
    subset_idx = cluster_idx ~= 0; % remove unassigned cells
    X_subset = X_macrophage( subset_idx, : );
    column_data_subset = column_data_macrophage( subset_idx,:);
    tsne_subset = tsne_macrophage(subset_idx,:);
    cluster_idx = cluster_idx( subset_idx );
   
    
    %% Compute AUC values
    
    cluster_names = unique(cluster_idx);
    num_clusters = length(cluster_names);
    num_genes = size(X_subset,2);

    AUC = zeros(num_genes,num_clusters);
    tic
    parfor num_cluster = 1:num_clusters
        cluster_ind = cluster_idx == cluster_names(num_cluster);
        for num_gene = 1:num_genes
            if ~mod(num_gene,100)
                fprintf('Cluster %i Gene %i of %i\n',num_cluster,num_gene,num_genes)
            end
            [~,~,~,AUC(num_gene,num_cluster)] = perfcurve( cluster_ind, X_subset(:,num_gene),true);
        end
    end
    toc
    
    %% Sort to find top markers
    [sort_AUC,sort_idx] = sort(AUC,'descend');    
    mkdir('results/macrophage_analysis/')
    for num_cluster = 1:num_clusters
        writetable(cell2table({row_data.GeneSymbol{sort_idx(1:100,num_cluster)}}','VariableNames', {strcat('cluster',num2str(num_cluster))}),...
                   strcat('results/macrophage_analysis/cluster',num2str(num_cluster),'.csv'), 'WriteVariableNames',false ); 

%         fprintf('Cluster %i\n',num_cluster);
%         for i = 1:10
%             fprintf('\t%s\n',row_data.GeneSymbol{sort_idx(i,num_cluster)})
%         end
    end
    writetable(cell2table(row_data.GeneSymbol,'VariableNames', {'GeneSymbol'}),...
                   strcat('results/macrophage_analysis/gene_symbols.csv'), 'WriteVariableNames',false ); 

    
    
end