function [ RL_list ] = import_RL_list( RL_path )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    %RL_path = '../data/receptor_ligand_Ramilowski.csv';
    RL_list = readtable( RL_path );
   % Use only literature supported interactions
    RL_list = RL_list( strcmp( RL_list.Pair_Evidence, 'literature supported') ,...
                       {'Pair_Name', 'Ligand_ApprovedSymobl', 'Receptor_ApprovedSymbol'});

   
end

