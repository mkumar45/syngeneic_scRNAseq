RL_pairs = interaction_scores{1}.Properties.RowNames;
split_names = cellfun( @(x) strsplit(x, '_'), RL_pairs, 'UniformOutput', false );
ligand_symbol = cellfun( @(x) x{1}, split_names, 'UniformOutput', false );
receptor_symbol = cellfun( @(x) x{2}, split_names, 'UniformOutput', false );


ligand_correlation = readtable( 'results/tirosh_phenotype/ligand_correlation_pct_treg_tfac.csv','ReadRowNames', true  ); 
receptor_correlation = readtable( 'results/tirosh_phenotype/receptor_correlation_pct_treg_tfac.csv','ReadRowNames', true  );

num_samples = length(ligand_expression);
num_cell_types = size(ligand_expression{1},2);

receptor_mat = zeros( [ size(ligand_expression{1}), num_samples] );
ligand_mat = zeros( size(receptor_mat) );


for n = 1:num_samples
    ligand_mat(:,:,n) = ligand_expression{n}{:,:};
    receptor_mat(:,:,n) = receptor_expression{n}{:,:};
end

[ unique_ligands, ligand_idx ] = unique( ligand_symbol, 'first' );
[ unique_receptors, receptor_idx ] = unique( receptor_symbol, 'first' );

ligand_subset = ligand_mat(ligand_idx,:,:);
receptor_subset = receptor_mat(receptor_idx,:,:);


%% Calculate all possible ligand/receptor scores (i.e. random pairs ) 
% num_unique_ligands = length(unique_ligands);
% num_unique_receptors = length(unique_receptors);
% permuted_scores = cell(num_cell_types,1);  
% for ct = 1:num_cell_types
%     permuted_scores{ct} = zeros( num_unique_ligands, num_unique_receptors, num_samples);
%     for l = 1:num_unique_ligands
%         permuted_scores{ct}(l,:,:) = bsxfun( @times, squeeze(ligand_subset(l,ct,:))', squeeze(receptor_subset(:,ct,:)) );
%     end
% end
% rho_correlations = cell(num_cell_types,1);  
% pVal_correlations = cell(num_cell_types,1);
% for ct = 1:num_cell_types
%     mat = reshape( permuted_scores{ct}, [ num_unique_ligands*num_unique_receptors,num_samples] );
%     [C, p_values] = corr( mat', phenotype_table{:,2} , 'type', 'spearman', 'rows', 'pairwise' );
%     rho_correlations{ct} = reshape( C, [ num_unique_ligands,num_unique_receptors] );
%     pVal_correlations{ct} = reshape( p_values, [ num_unique_ligands,num_unique_receptors] );
% end 

%% Calculate all possible ligand/receptor scores (i.e. random pairs ) 
num_unique_ligands = length(unique_ligands);
num_unique_receptors = length(unique_receptors);
full_permuted = zeros(  num_cell_types*num_unique_ligands, num_cell_types*num_unique_receptors, num_samples );
for ct_ligand = 1:num_cell_types
    ct_ligand
    for ct_receptor = 1:num_cell_types
        permuted_scores = zeros( num_unique_ligands, num_unique_receptors, num_samples);
        for l = 1:num_unique_ligands
            permuted_scores(l,:,:) = bsxfun( @times, squeeze(ligand_subset(l,ct_ligand,:))', squeeze(receptor_subset(:,ct_receptor,:)) );
        end
        full_permuted( (ct_ligand-1)*num_unique_ligands+1: ct_ligand*num_unique_ligands,...
                       (ct_receptor-1)*num_unique_receptors+1: ct_receptor*num_unique_receptors, : ) = permuted_scores;             
    end 
end

tic
mat = reshape( full_permuted, [ num_cell_types^2*num_unique_ligands*num_unique_receptors,num_samples] );
[rho_correlations, pVal_correlations] = corr( mat', phenotype_table{:,2} , 'type', 'spearman', 'rows', 'pairwise' );
toc

%% Plot
ligand_correlation_subset = ligand_correlation(ligand_idx,:);
receptor_correlation_subset = receptor_correlation(receptor_idx,:);

bins = -1:.1:1;
ct=5;
var = 1
for i = 1:length(bins)-1
    figure(1)
    subplot(10,2,i)
    idx = ligand_correlation_subset{:,ct} < bins( i + 1) & ligand_correlation_subset{:,ct} > bins(i);
    corr_vec = rho_correlations{ct}(idx,:);    
    
    is_gt = bsxfun( @gt, abs(corr_vec), abs(ligand_correlation_subset{idx,ct}) );
    
    histogram( corr_vec(:),40 )
    med = median( ligand_correlation_subset{idx,ct} );
    line([med, med], ylim, 'LineWidth', 2, 'Color', 'r');
    title( sprintf( 'Bin %1.1f to %1.1f', bins(i), bins(i+1) ) )
    xlim([-1,1])
    if ~isempty( corr_vec(:) )
        [f,xi] = ksdensity(corr_vec(:));
        figure(2)
        subplot(6,2,var)
        plot(xi,f,'k','Linewidth',2)
        line([med, med], ylim, 'LineWidth', 2, 'Color', 'r');
        xlim([-1,1])
        xlabel('Correlation with random receptors'); ylabel('Density')
        title( sprintf( 'Ligand correlation between %1.1f to %1.1f', bins(i), bins(i+1) ) )
        text(-0.8, 1, sprintf( ' p < %0.3f ', sum(sum(is_gt))./numel(is_gt) ) )
        var = var + 1;
    end
end
    
%% Receptor
ct=3;
var = 1
for i = 1:length(bins)-1
    figure(3)
    subplot(10,2,i)
    idx = receptor_correlation_subset{:,ct} < bins( i + 1) & receptor_correlation_subset{:,ct} > bins(i);
    corr_vec = rho_correlations{ct}(:,idx);   
    

    
    histogram( corr_vec(:),40 )
    med = median( receptor_correlation_subset{idx,ct} );
    line([med, med], ylim, 'LineWidth', 2, 'Color', 'r');
    title( sprintf( 'Bin %1.1f to %1.1f', bins(i), bins(i+1) ) )
    xlim([-1,1])
    if ~isempty( corr_vec(:) )
        [f,xi] = ksdensity(corr_vec(:));
        figure(4)
        subplot(7,2,var)
        plot(xi,f,'k','Linewidth',2)
        line([med, med], ylim, 'LineWidth', 2, 'Color', 'r');
        xlim([-1,1])
        is_gt = bsxfun( @gt, abs(corr_vec), abs(receptor_correlation_subset{idx,ct})' );
        xlabel('Correlation with random ligand'); ylabel('Density')
        title( sprintf( 'Receptor correlation between %1.1f to %1.1f', bins(i), bins(i+1) ) )
        text(-0.8, 1, sprintf( ' p < %0.3f ', sum(sum(is_gt))./numel(is_gt) ) )
        var = var + 1;
    end
end

