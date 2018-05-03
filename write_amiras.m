addpath('C:/Users/DTHUFF/Documents/Research/generalCodes/');

d = '//mpufs5/data_wnx1/_Data/FDG Bone/';

patients = {'01','02','03','04','05','06','07','08','09','10',...
            '11','12','13','14','15','16','17','18','19','20'};

%% HB: CT + PET Recon 1
for i=1:size(patients,2)
    
    d_dicom_ct = [d 'HB' patients{i} '/Raw/CT/'];
    %d_dicom_pet = [d 'HB' patients{i} '/Raw/Recon1/'];
    
    d_amira_ct = [d 'HB' patients{i} '/Processed/HB' patients{i} '_CT_no_table.am'];
    %d_amira_pet = [d 'HB' patients{i} '/Processed/HB' patients{i} '_FDG_SUV.am'];
    
    ct_struct = dcm2mat(d_dicom_ct);
    ct_struct.data = removeTable(ct_struct.data);
    mat2am(ct_struct, d_amira_ct);
    
    %pet_struct = dcm2mat(d_dicom_pet);
    %mat2am(pet_struct, d_amira_pet);
end
    
    
%% AB: CT + PET Recon 1
for i=1:10
    
    d_dicom_ct = [d 'AB' patients{i} '/Raw/CT/'];
    d_dicom_pet = [d 'AB' patients{i} '/Raw/Recon1/'];
    
    d_amira_ct = [d 'AB' patients{i} '/Processed/AB' patients{i} '_CT.am'];
    d_amira_pet = [d 'AB' patients{i} '/Processed/AB' patients{i} '_FDG_SUV.am'];
    
    ct_struct = dcm2mat(d_dicom_ct);
    mat2am(ct_struct, d_amira_ct);
    
    pet_struct = dcm2mat(d_dicom_pet);
    mat2am(pet_struct, d_amira_pet);
end




