function [ output_args ] = plot_slice( scan, slice )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    figure();
    if max(scan(:)) == 1 %for masks
        imshow(squeeze(scan(:,:,slice)));
    else
        imshow(uint8(squeeze(scan(:,:,slice))));
    end        
end

