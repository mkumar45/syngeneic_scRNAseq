function [ cMap_ct, cMap_model, red_blue_cmap ] = color_maps()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    cMap_ct = [ 166,206,227 % B16F10 1
                177,89,40, % Bcell 2
                31,120,180 % CAF 3
                178,223,138 % CT-26 4
                51,160,44 % EMT6 5 
                0,0,0 % Endothelial 6 
                251,154,153 % LL2 7 
                227,26,28 % Macrophage 8
                253,191,111 % MC38 9
                255,127,0 % NK 10 
                202,178,214 % SA1N 11 
                106,61,154 % Tcell 12 
                125,125,125% unassigned 13 
                255,255,153]./256; %pDC 14
            
            
    cMap_model = [  166,206,227 % B16F10
                    178,223,138 % CT-26
                    51,160,44 % EMT6
                    251,154,153 % LL2
                    253,191,111 % MC38
                    202,178,214 ]./256; % SA1N
                
                
    color_vec = transpose(0:.02:1);
    zero_vec = zeros(size(color_vec));
    red_blue_cmap = [ zero_vec zero_vec color_vec(end:-1:1); color_vec zero_vec zero_vec ];


end

