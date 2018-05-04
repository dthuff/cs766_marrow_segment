%% plotting
    % voxel aspect ratio, for plotting
    voxel_ratio = ct_struct.voxel_size(3)/ct_struct.voxel_size(2);

    % generic axial space vector
    z_ = 1:size(ct, 3);

    %% Sagittal images showing bone masking workflow
    fh1 = figure('position',[50, 200, 1800, 800]);
    [~, slice] = max(sum(sum(whole_bone_mask,1), 3)); 
    
    %sagittal slice with most bone voxels
    slice_to_plot = squeeze(ct(:,slice,:));
    ct_window = [-700 800];
    pet_window = [0 8];

    % Raw CT
    subplot(1, 4, 1);
    imshow(rot90(slice_to_plot), ct_window);
    daspect([voxel_ratio 1 1]);
    title(['HB' patients{i} ' Raw CT Slice']);

    % After HU Threshold
    subplot(1, 4, 2);
    imshow(rot90(slice_to_plot > bone_threshold));
    daspect([voxel_ratio 1 1]);
    title(['Bone Threshold HU > ' num2str(bone_threshold)]);

    % After Mask Cleanup
    subplot(1, 4, 3);
    imshow(rot90(squeeze(whole_bone_mask(:,slice,:))));
    daspect([voxel_ratio 1 1]);
    title('Whole Bone Mask After Preprocessing');

    % Masked CT Annotated with found minima
    subplot(1, 4, 4);
    ct_bone(~ct_bone) = -1000;
    %imshow(rot90(squeeze(ct_bone(:,slice,:))), ct_window);
    imshow(rot90(squeeze(pet_bone(:,slice,:))), pet_window);

    hold on;
    daspect([voxel_ratio 1 1]);
    title('Raw CT*Whole Bone Mask');
    other_marrow_mask = saveAnnotatedImg(fh1);
    imwrite(other_marrow_mask, [d_proj 'images/HB' patients{i} '_bone_mask.png']);


    %% Axial mean HU, SUV with minima locations
    fh2 = figure('Position',[50, 500, 1000, 500]);
    % CT
    subplot(2, 1, 1);
    plot(z_, mean_hu);
    hold on;
    plot(ct_minima, mean_hu(ct_minima), '*');
    ylim([100,500]);
    xlim([shoulder_slice, pelvis_slice]);
    title(['HB' patients{i} ' Axial Mean Image Intensity Values']);
    ylabel('Mean HU');
    xlabel('z');

    % PET
    subplot(2, 1, 2);
    plot(z_, mean_suv);
    hold on;
    plot(pet_minima, mean_suv(pet_minima), '*');
    ylim([0, 3]);
    xlim([shoulder_slice, pelvis_slice]);
    ylabel('Mean SUV');
    xlabel('z');
    other_marrow_mask = saveAnnotatedImg(fh2);
    %imwrite(other_marrow_mask, [d_proj 'images/HB' patients{i} '_axial_dists.png']);


    %% CT, PET sagittal slices annotated with identified minima locations
    fh3 = figure('position', [50, 300, 900, 600]);
    subplot(1, 2, 1);
    imshow(rot90(slice_to_plot), ct_window);
    hold on;
    %visboundaries(rot90(squeeze(bone_marrow_mask(:,slice,:))), 'linewidth',1, 'color','blue');
    for j=1:size(ct_minima_inv,2)
        line([300, 310], [ct_minima_inv(j), ct_minima_inv(j)], 'color', 'red', 'linewidth',1);
    end
    for j=1:size(pet_minima_inv,2)
        line([310, 320], [pet_minima_inv(j), pet_minima_inv(j)], 'color', 'green', 'linewidth',1);
    end
    daspect([voxel_ratio 1 1]);
    title(['HB' patients{i} ' CT Slice with Detected Minima']);
    
    subplot(1, 2, 2);
    imshow(rot90(squeeze(pet(:,slice,:))), pet_window);
    hold on;
    %visboundaries(rot90(squeeze(bone_marrow_mask(:,slice,:))), 'linewidth',1, 'color','blue');
    for j=1:size(ct_minima_inv,2)
        line([300, 310], [ct_minima_inv(j), ct_minima_inv(j)], 'color', 'red', 'linewidth',1);
    end
    for j=1:size(pet_minima_inv,2)
        line([310, 320], [pet_minima_inv(j), pet_minima_inv(j)], 'color', 'green', 'linewidth',1);
    end
    daspect([voxel_ratio 1 1]);
    title(['HB' patients{i} ' PET Slice with Detected Minima']);
    
    other_marrow_mask = saveAnnotatedImg(fh3);
    %imwrite(other_marrow_mask, [d_proj 'images/HB' patients{i} 'marrow2.png']);

    
    %% whole body mask
    figure()
    imshow(rot90(squeeze(body_mask(:, 256, :))))
    daspect([voxel_ratio 1 1]);


    %% shoulder, pelvis localization workflow
    figure('Position',[200,300,800,600]);
    subplot(3, 1, 1);
    plot(z_, body_mask_axial);
    title('Body X-Sectional Area A');

    subplot(3, 1, 2);
    plot(z_(1:end-1), abs(diff(body_mask_axial)));
    title('A'' = |dA/dz|');

    subplot(3, 1, 3);
    plot(z_(1:end-1), d_body_mask_axial);
    hold on;
    title('A'' \ast N(0, \sigma^2)');
    xlabel('z');
    ylim([0,12000]);
    plot([pelvis_slice, pelvis_slice],[0,12000],'k--');
    plot([shoulder_slice-shoulders_offset, shoulder_slice-shoulders_offset],[0,12000],'k--');

    
    %% other marrow workflow
    figure('position',[50, 200, 1800, 800]);

    % vertebral marrow mask
    subplot(1, 4, 1);
    imshow(rot90(squeeze(vertebra_marrow_mask(:,slice,:))));
    daspect([voxel_ratio 1 1]);
    title('Final Vertebral Marrow');
    
    % bone mask - vertebral marrow mask
    subplot(1, 4, 2);
    imshow(rot90(squeeze(other_bone_mask(:,slice,:))));
    daspect([voxel_ratio 1 1]);
    title('Bone Mask - Vertebral Marrow');
    
    % after erosion
    subplot(1, 4, 3);
    imshow(rot90(squeeze(bone_mask_eroded(:,slice,:))));
    daspect([voxel_ratio 1 1]);
    title('Eroded Difference Mask');
    
    % largest components
    subplot(1, 4, 4);
    imshow(rot90(squeeze(other_marrow_mask(:,slice,:))));
    daspect([voxel_ratio 1 1]);
    title('Largest Other Marrow Components');
    