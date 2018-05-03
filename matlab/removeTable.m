function CT_tableRemoved = removeTable(CT_original)

CT_resamp = imresize3(CT_original, size(CT_original)./3);

% Find areas of air and perform a closing and opening
air = CT_resamp<-700;
se = strel('sphere', 3);
air = imclose(air, se); 
air = imopen(air, se);

% Find largest component that is not air - this will be the patient
CC = bwconncomp(~air);
sizeCC = zeros([1, CC.NumObjects]);
for conc = 1:CC.NumObjects
    sizeCC(conc) = size(CC.PixelIdxList{conc},1);
end
[~,I] = sort(sizeCC, 'descend');
concCList = CC.PixelIdxList(:,I);
patient_mask = zeros(size(CT_resamp));

% Close holes and resize back to original CT size
patient_mask(concCList{1}) = 1;
patient_mask = imfill(patient_mask, 'holes');
se = strel('sphere', 6);
patient_mask = imclose(patient_mask, se);
patient_mask = imresize3(patient_mask, size(CT_original), 'nearest');

CT_tableRemoved = CT_original; 
CT_tableRemoved(patient_mask==0)=-1000;