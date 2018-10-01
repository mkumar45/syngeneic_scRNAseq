function [ return_symbols ] = find_homologs(search_symbols, search_species )
% Find human-mouse homologs
%[ return_ID ] = find_homologs(search_symbols, search_species )
%   INPUTS:
%       search_symbols - gene symbols to be converted. 
%       search_species - species of input symbols. Either 'human' or
%       'mouse, laboratory'
%   OUTPUTS: 
%       return_symbols: symbols of identified homologs


    homolog_table = readtable('data/HOM_MouseHumanSequence.rpt.txt',...
                              'ReadVariableNames',true,...
                              'ReadRowNames',false,...
                              'Delimiter','\t');

    search_organism_idx = strcmp(homolog_table{:,'CommonOrganismName'},search_species);
    find_input = @(input_symbol) ismember(homolog_table{:,'Symbol'},input_symbol) & search_organism_idx;
    get_homolog_ID = @(input_idx) ismember(homolog_table{:,'HomoloGeneID'}, homolog_table{input_idx,'HomoloGeneID'});
    get_homolog_idx = @(homolog_idx)  homolog_idx & ~search_organism_idx ;
    get_homologs = @(return_idx) homolog_table{return_idx,'Symbol'} ;

    return_symbols = cellfun( @(x) get_homologs(get_homolog_idx(get_homolog_ID(find_input(x)))), search_symbols, 'UniformOutput', false );
    
end

