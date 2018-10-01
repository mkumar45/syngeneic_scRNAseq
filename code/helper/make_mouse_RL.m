function make_mouse_RL( RL_path )
%MAKE_MOUSE_RL Summary of this function goes here
%   Detailed explanation goes here
    %RL_path = '../';
    RL_list = readtable( RL_path );
   % Use only literature supported interactions
    RL_list = RL_list( strcmp( RL_list.Pair_Evidence, 'literature supported') , :);

                   
   mouse_ligands = find_homologs( RL_list.Ligand_ApprovedSymbol, 'human' );
   mouse_receptors = find_homologs( RL_list.Receptor_ApprovedSymbol, 'human' );
   
   unique_func = @(x) length(x) == 1;
   unique_ligand = cellfun( unique_func, mouse_ligands);
   unique_receptor = cellfun( unique_func, mouse_receptors);
   
   interaction_idx = unique_ligand & unique_receptor;
   
   mouse_RL_list = RL_list(interaction_idx,:);
   
   unravel_func = @(x) x{1};
   mouse_RL_list.Ligand_ApprovedSymbol = cellfun( unravel_func, mouse_ligands(interaction_idx), 'UniformOutput',false );
   mouse_RL_list.Receptor_ApprovedSymbol = cellfun( unravel_func, mouse_receptors(interaction_idx), 'UniformOutput',false );
   
   writetable(mouse_RL_list,'data/mouse_receptor_ligand.csv')
   
end

