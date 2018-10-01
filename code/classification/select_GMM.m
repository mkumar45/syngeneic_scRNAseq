function [ best_models, cluster_idx, sort_idx, parameter_tbl ] = select_GMM( X, models, gene_names )
%SELECT_GMM Select number of components for each GMM model 
%   INPUTS:
%       X: expression data (used for clustering after GMM is selected
%       models: cell containing GMMs. Each row should match genes in X.  
%   OUTPUTS:
%       best_models: GMM for each gene fit to all data using selected
%       number of components
%       cluster_idx: cluster labels for each cell using GMM
%       sort_idx: sorted order of mixture componenets in descending order
%       by mean

    options = statset('MaxIter',500);
    % Extract AIC and BIC from fit mixture models
    AIC = cellfun( @(x) x.AIC , models );
    BIC = cellfun( @(x) x.BIC , models );
    % Compute mean + standard deviation across k-fold partitions
    mean_BIC = mean(BIC,3) ;
    std_BIC = std(BIC, [], 3);

    num_genes = size( models,1);
    best_models = cell(num_genes,1);
    cluster_idx = zeros( size(X,1), num_genes );
    sort_idx = cell( 1, num_genes);
    
    % Select number of conponents within 3x standard error of the mean of
    % the minimum value
    for n = 1:num_genes
        [~,min_idx] = min(mean_BIC(n,:));
        num_comp = find( mean_BIC(n,:) - 3./sqrt(5).*std_BIC(n,:) <  mean_BIC(n,min_idx) + 3./sqrt(5)*std_BIC(n,min_idx) , 1, 'first' ) ;

        best_models{n} = fitgmdist( X( : , n ), num_comp, 'Options',options,...
                                    'CovarianceType','diagonal','RegularizationValue',0.01);
        cluster_idx( :, n ) = cluster(best_models{n}, X( :, n) );
        
        %Plot GMM
        figure
        subplot(121)
        x_range = -1:.1:5;
        norm_pdf = normpdf( x_range, best_models{n}.mu, squeeze(best_models{n}.Sigma ) );
        for c = 1:num_comp
            plot( x_range, norm_pdf(c,:))
            hold on
        end
        axis square; box on;
        subplot(122)
        errorbar( 1:5, mean_BIC(n,:), 3./sqrt(5).*std_BIC(n,:) ,'k' )
        hold on
        plot( num_comp, mean_BIC(n,num_comp), 'r.')
        xlim([0,6]);
        xlabel('Number of mixture components'); ylabel('Bayesian Info Criteria');
        set(gca,'FontSize',16)
        set(gcf, 'PaperPositionMode', 'auto'); legend('boxoff','Location','NorthEast') 
        axis square; box on;
        suptitle(gene_names(n))
        
        [~,sort_idx{n}] = sort(best_models{n}.mu,'descend');
        
    end
    
    num_comp = cell2table(cellfun( @(x) length(x.mu), best_models, 'UniformOutput', false ),'VariableNames',{'Num_components'});
    mu_tbl = cell2table(cellfun( @(x) x.mu, best_models, 'UniformOutput', false ),'VariableNames',{'Mu'});
    sigma_tbl = cell2table(cellfun( @(x) squeeze(x.Sigma), best_models, 'UniformOutput', false ),'VariableNames',{'Sigma'});

    parameter_tbl = [num_comp, mu_tbl, sigma_tbl];
    parameter_tbl.Properties.RowNames =  gene_names;
    
end

