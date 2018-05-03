function [shoulder_slice, pelvis_slice] = find_spine_ends(pet, shoulder_offset, pelvis_offset)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    body_mask = pet > 0.1;
    
    % sum over each axial slice
    body_mask_axial = squeeze(sum(sum(body_mask, 1), 2));
    
    % gaussian smoothed, abs, 1st derivative
    w = gausswin(20);
    d_body_mask_axial = filter(w, 1, abs(diff(body_mask_axial)));
    
    % biggest peaks should be ~shoulders and ~pelvis due to change in x
    % sectional area of body
    [~, locs, ~, p]  = findpeaks(d_body_mask_axial);
    [~, ind] = sort(p, 'descend');
    
    pelvis_slice = max(locs(ind(1:2)))+pelvis_offset;
    shoulder_slice = min(locs(ind(1:2)))+shoulder_offset;
end
