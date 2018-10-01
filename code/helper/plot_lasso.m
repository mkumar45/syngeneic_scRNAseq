function [ output_args ] = plot_lasso( B, stats, fName )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    figure('units','normalized','outerposition',[0 0 1 1])

    lassoPlot(B,stats,'PlotType','CV');
    box on; axis square; set(gca,'FontSize',16)
    
    min_var_names = stats.PredictorNames( B( :, stats.IndexMinMSE ) ~= 0 )'
    se_var_names = stats.PredictorNames( B( :, stats.Index1SE ) ~= 0 )'
    
    figure
    subplot(121)    
    bar( 1:length(min_var_names),B( B( :, stats.IndexMinMSE ) ~= 0 ) )
    set(gca,'TickLabelInterpreter','none');
    set(gca,'XTickLabels', min_var_names); set(gca,'XTickLabelRotation',45);
    box on; axis square; set(gca,'FontSize',10)

    subplot(122)
    bar( 1:length(se_var_names), B( B( :, stats.Index1SE ) ~= 0 ) )
    set(gca,'TickLabelInterpreter','none');
    set(gca,'XTickLabels', se_var_names); set(gca,'XTickLabelRotation',45);
    box on; axis square; set(gca,'FontSize',10)
	%saveas(gcf, strcat('figures/',fName,'.pdf'));

end