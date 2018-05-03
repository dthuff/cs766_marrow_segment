clear all;

addpath('C:/Users/DTHUFF/Documents/Research/generalCodes/');

d_proj = 'C:/Users/DTHUFF/Documents/Class Material/CS766/project/';

patients = {'01','02','03','04','05','06','07','08','09','10',...
            '11','12','13','14','15','16','17','18','19','20'};

for i=1:size(patients, 2)
    tic;
%% Load Patient Data

    d_ct = ['//mpufs5/data_wnx1/_Data/FDG Bone/HB' patients{i} '/Raw/CT'];
    d_pet = ['//mpufs5/data_wnx1/_Data/FDG Bone/HB' patients{i} '/Raw/Recon1'];

    [ct_struct, ct_header] = dcm2mat(d_ct);
    [pet_struct, pet_header] = dcm2mat(d_pet);

%% preprocessing
    ct = ct_struct.data;
    pet = imresize3(pet_struct.data, size(ct));

    % extract bone mask
    bone_threshold = 150;
    small_volume_threshold = 1000;
    dilation_scale = 3;

    whole_bone_mask = extract_bone_mask(ct, bone_threshold, small_volume_threshold, dilation_scale);

    % mask CT, PET with bone mask
    ct_bone = ct.*whole_bone_mask;
    pet_bone = pet.*whole_bone_mask;
    vertebra_marrow_mask = zeros(size(ct));

    %% segment vertebral marrow 
    % find axial minima on both CT and PET
    ct_prominence_threshold = 20;
    pet_prominence_threshold = 0.15;

    m_ct = find_axial_minima(ct, whole_bone_mask, ct_prominence_threshold);
    ct_minima = fliplr(m_ct.minima_slice_indices);
    mean_hu = m_ct.means;

    m_pet = find_axial_minima(pet, whole_bone_mask, pet_prominence_threshold);
    pet_minima = fliplr(m_pet.minima_slice_indices);
    mean_suv = m_pet.means;    
    
    % find pelvis, shoulders, only look for vertebra between those limits
    shoulder_offset = 50;
    pelvis_offset = 0;
    
    [shoulder_slice, pelvis_slice] = find_spine_ends(pet, shoulder_offset, pelvis_offset);
    
    pet_minima = pet_minima((pet_minima < pelvis_slice) & (pet_minima > shoulder_slice));
    ct_minima =  ct_minima((ct_minima < pelvis_slice) & (ct_minima > shoulder_slice));

    % inverted minima indices, for plotting
    pet_minima_inv = size(ct,3) + 1 - pet_minima;
    ct_minima_inv =  size(ct,3) + 1 - ct_minima;

    
    % Kalman filtering?
    
    
    
    % for each pet minima
    for j = 1:size(pet_minima, 2)-1
        
        %crop 2 pixels off top and bottom
        vertebra_start = pet_minima(j+1) + 2;
        vertebra_stop = pet_minima(j) - 2;
        
        this_vertebra_mask = m_ct.mask(:,:,vertebra_start:vertebra_stop);
        
        % erosion
        this_vertebra_mask = imerode(this_vertebra_mask, strel('sphere', 3));
        
        % grab largest blob if more than one component is present
        if max(max(max(bwlabeln(this_vertebra_mask)))) > 1
            this_vertebra_mask = ExtractNLargestBlobs3(this_vertebra_mask, 1);
        end   
        %this_vertebra_mask = activecontour(pet(:,:,start:stop), this_vertebra_mask, 100);
        
        vertebra_marrow_mask(:,:,vertebra_start:vertebra_stop) = vertebra_marrow_mask(:,:,vertebra_start:vertebra_stop) + this_vertebra_mask;
    end
    
    other_bone_mask = whole_bone_mask-vertebra_marrow_mask;
    
    bone_mask_eroded = imerode(other_bone_mask, strel('sphere',4));
    
    n_marrow_blobs = 5;
    other_marrow_mask = bwlabeln(ExtractNLargestBlobs3(bone_mask_eroded, n_marrow_blobs));
    
    bone_marrow_mask = (vertebra_marrow_mask > 0) | (other_marrow_mask > 0);
    
%% plotting
    % Ground truth axial slices for patient HB03
    % gt_slices = 336-[95, 104, 111, 118, 126, 135, 144, 153, 164, 174, 185, 196];

    % voxel aspect ratio, for plotting
    voxel_ratio = ct_struct.voxel_size(3)/ct_struct.voxel_size(2);

    % generic axial space vector
    z_ = 1:size(ct, 3);

    %% Sagittal images showing bone masking workflow
    fh1 = figure('position',[50, 200, 1800, 800]);
    [~, slice] = max(sum(sum(whole_bone_mask,1), 3)); %sagittal slice with most bone voxels
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
    imwrite(other_marrow_mask, ['./images/HB' patients{i} '_bone_mask.png']);


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
    imwrite(other_marrow_mask, ['./images/HB' patients{i} '_axial_dists.png']);


    %% CT, PET sagittal slices annotated with identified minima locations
    fh3 = figure('position', [50, 300, 900, 600]);
    subplot(1, 2, 1);
    imshow(rot90(slice_to_plot), ct_window);
    hold on;
    %visboundaries(rot90(squeeze(bone_marrow_mask(:,slice,:))), 'linewidth',1, 'color','blue');
% 
%     for j=1:size(ct_minima_inv,2)
%         line([300, 310], [ct_minima_inv(j), ct_minima_inv(j)], 'color', 'red', 'linewidth',1);
%     end
%     for j=1:size(pet_minima_inv,2)
%         line([310, 320], [pet_minima_inv(j), pet_minima_inv(j)], 'color', 'green', 'linewidth',1);
%     end
    daspect([voxel_ratio 1 1]);
    title(['HB' patients{i} ' CT Slice with Detected Minima']);
    
    subplot(1, 2, 2);
    imshow(rot90(squeeze(pet(:,slice,:))), pet_window);
    hold on;
    %visboundaries(rot90(squeeze(bone_marrow_mask(:,slice,:))), 'linewidth',1, 'color','blue');
%     for j=1:size(ct_minima_inv,2)
%         line([300, 310], [ct_minima_inv(j), ct_minima_inv(j)], 'color', 'red', 'linewidth',1);
%     end
%     for j=1:size(pet_minima_inv,2)
%         line([310, 320], [pet_minima_inv(j), pet_minima_inv(j)], 'color', 'green', 'linewidth',1);
%     end
    daspect([voxel_ratio 1 1]);
    title(['HB' patients{i} ' PET Slice with Detected Minima']);
    
    other_marrow_mask = saveAnnotatedImg(fh3);
    imwrite(other_marrow_mask, ['./images/HB' patients{i} 'marrow2.png']);

    
    %%
    
    figure()
    imshow(rot90(squeeze(body_mask(:, 256, :))))
    daspect([voxel_ratio 1 1]);


    %%
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

    
    %%
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
    
    %six largest components
    subplot(1, 4, 4);
    imshow(rot90(squeeze(other_marrow_mask(:,slice,:))));
    daspect([voxel_ratio 1 1]);
    title('Largest Other Marrow Components');
    
    
%% write out to amira
    bm_out.data = single(other_marrow_mask);
    bm_out.voxel_size = ct_struct.voxel_size;
    bm_out.start = ct_struct.start;
    mat2am(bm_out, ['//mpufs5/data_wnx1/_Data/FDG Bone/HB' patients{i} '/Processed/HB' patients{i} '_other_eroded_largest_BM_MASK.am']);
    
%%  Close figures, write out timing
    close all;
    
    t=toc;
    disp(['Done with patient HB' patients{i} ' in ' num2str(round(t)) ' seconds...']);    
end





