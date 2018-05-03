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

    whole_bone_mask = extract_bone_mask(ct, bone_threshold, ...
                            small_volume_threshold, dilation_scale);

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
    
    pet_minima = pet_minima((pet_minima < pelvis_slice) & ...
                            (pet_minima > shoulder_slice));
    ct_minima =  ct_minima((ct_minima < pelvis_slice) & ...
                           (ct_minima > shoulder_slice));

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
        
        vertebra_marrow_mask(:,:,vertebra_start:vertebra_stop) = ...
            vertebra_marrow_mask(:,:,vertebra_start:vertebra_stop) + this_vertebra_mask;
    end
    
    % subtract vertebra, erode cortical bone
    other_bone_mask = whole_bone_mask-vertebra_marrow_mask;
    
    bone_mask_eroded = imerode(other_bone_mask, strel('sphere',4));
    
    % grab largest components as other marrow of interest 
    n_marrow_blobs = 5;
    other_marrow_mask = bwlabeln(ExtractNLargestBlobs3(bone_mask_eroded, n_marrow_blobs));
    
    %final seg = union of vertebral, other marrow
    bone_marrow_mask = (vertebra_marrow_mask > 0) | (other_marrow_mask > 0);
    
%% write out to amira
    bm_out.data = single(bone_marrow_mask);
    bm_out.voxel_size = ct_struct.voxel_size;
    bm_out.start = ct_struct.start;
    mat2am(bm_out, ['//mpufs5/data_wnx1/_Data/FDG Bone/HB' patients{i} '/Processed/HB' patients{i} '_BM_MASK.am']);
    
%%  Close figures, write out timing
    close all;
    
    t=toc;
    disp(['Done with patient HB' patients{i} ' in ' num2str(round(t)) ' seconds...']);    
end





