function [ indices ] = find_indices( cell_array_1, cell_array_2 )
%FIND_INDICES Find the indices of elements cell_array_1 in cell_array_2
%   INPUTS:
%       cell_array_1 - cell array containing items you want to search for
%       in cell_array_2
%       cell_array_2 - 
%   OUTPUTS:
%       indices - array of length(cell_array_1). indices(i) contains the
%       index of cell_array_1(i) in cell_array_2
    
    unique_vals = unique(cell_array_1);
    indices = zeros(size(cell_array_1,1),1);
    for val = 1:length(unique_vals)
       loc = find(strcmp(cell_array_2,unique_vals{val}));  % duplicate genes (fix) 
       if ~isempty(loc)
           indices( strcmp( cell_array_1, unique_vals{val} ) ) = loc(1); 
       else
           unique_vals{val};
       end
    end
    
end
