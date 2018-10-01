function [ marker_table ] = read_marker_table( marker_path )
%READ_MARKER_TABLE Import a marker table specifiying cell-types and gene markers 
%   INPUTS:
%       marker_path: filepath to user definied marker path  
%   OUTPUTS:
%       marker_table: matlab table containing cell-type and marker
%       information
    marker_table = readtable( marker_path, 'ReadRowNames',true, 'ReadVariableNames',true);   
    marker_table{:,:} = regexprep(marker_table{:,:},'AND','1'); 
    marker_table{:,:} = regexprep(marker_table{:,:},'NOT','0');
    marker_table{:,:} = regexprep(marker_table{:,:},'OR','-1');
    marker_table.Properties.VariableNames = regexprep(marker_table.Properties.VariableNames,'_','-'); % gene names are converted when importing table

end

