function [ mask ] = extract_bone_mask( ct, bone_threshold, small_volume_threshold, dilation_scale )
    % extract_bone_mask 
    
    % remove table
    ct = removeTable(ct);
    
    % threshold at hu > bone_threshold;
    mask = ct > bone_threshold;

    % reject small volumes
    mask = bwareaopen(mask, small_volume_threshold);

    % fill holes
    mask = imfill(mask, 18, 'holes');

    % dilate and erode to fill crevices
    se = strel('sphere', dilation_scale);
    mask = imclose(mask, se); 
    mask = imopen(mask, se);
    
    
%     for i=1:size(ct, 1);
%         mask(i,:,:) = bwmorph(squeeze(mask(i,:,:)), 'dilate', dilation_scale);
%         mask(i,:,:) = bwmorph(squeeze(mask(i,:,:)), 'erode', dilation_scale);
%     end

    % fill holes again
    mask = imfill(mask, 18, 'holes');
end

