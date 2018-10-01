function [ output_args ] = write_null_scores( output_path, null_scores )
%WRITE_NULL_SCORES Summary of this function goes here
%   Detailed explanation goes here


    num_reps = size(null_scores,3);
    mkdir( output_path ) 

    for n = 1:num_reps
        
        csvwrite( strcat(output_path, 'null_',num2str(n),'.csv'),null_scores(:,:,3) );
    end


end