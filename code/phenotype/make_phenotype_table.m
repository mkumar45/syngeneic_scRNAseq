function make_phenotype_table( cell_type_path, column_path, growth_path, output_path )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



    [ tcell_table ] = compute_t_cell_pct( cell_type_path, column_path );
    growth_rates = import_growth_rate( growth_path );

    tcell_table.Model = tcell_table.Properties.RowNames;
    growth_rates.Model = growth_rates.Properties.RowNames;
    
    tbl = outerjoin(tcell_table, growth_rates, 'LeftVariables','pct_tcell',...
                                               'RightVariables','growth_rate');
    writetable( tbl, output_path );
end

