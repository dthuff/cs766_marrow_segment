function [ m ] = find_axial_minima( scan, mask, prominence_threshold )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    z_ = 1:size(scan, 3);

    means = zeros(size(z_));
    mins = zeros(size(z_));
    mask_out = zeros(size(scan));
    
    for i=1:size(z_, 2)
        
        this_slice = scan(:, :, i);
        this_mask_slice = mask(:, :, i);
        
        if any(any(this_mask_slice))
            
            this_mask_slice = ExtractNLargestBlobs(this_mask_slice, 1);
            
            mask_out(:,:,i) = this_mask_slice;
            means(i) = mean(this_slice(this_mask_slice)); 
            mins(i) = min(this_slice(this_mask_slice));
        end
    end

    [~, ind_min, ~, prominence] = findpeaks(-means);
    
    m.minima_slice_indices = ind_min(prominence>prominence_threshold);
    m.means = means;
    m.mins = mins;
    m.mask = mask_out;
end

