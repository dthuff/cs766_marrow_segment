function [Geometry, header] = dcm2mat(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Summary: 
% Reads in a directory of DICOM images and converts them into a Matlab
% image object. There must be no other files or folders in the directory. 
% It is similar to createFrames and dcm2am. Works for PET, CT, and MR. By
% default, the image is converted to SUV for PET, but includes an SUV 
% conversion factor.
% 
% Input:
%   directory = (string) directory containing the DICOM files. 
%               Eg, 'D:\PET\Image1' (optional)
% 
% Output:
%   Geometry = (structure) contains three (CT) or four (PET) fields
%   Geometry.data = (array) 3D image data in standard dicom unit (eg Bq/cc)
%   Geometry.voxel_size = (1x3 array) voxel sizes in mm for each dimension
%   Geometry.start = (1x3 array) physical starting location for each
%                    dimension in mm
%   Geometry.SUV_factor = (double) factor that converted Bq/cc into SUV. To
%                    convert back to Bq/cc, divide by this factor. Only for
%                    PET images -- not included in CT images.
%   header  = (structure) DICOM header of the first slice in the image
%                    series (optional)
% 
%Example usage:
% [G, H] = dcm2mat('C:\DICOM\Scan1');
% This read in the DICOMs in the Scan1 directory to an image object G and
% the header object H.
%
% Non-standard subroutines:
%   dicomread_dcm2mat.m 
%   Image processing toolbox
%   
%
%Tips:
% Very similar to dcm2am.m in function
% 
%Changes from previous versions:
% -Corrects for an earlier error in the time-decay correction
% -Doesn't exclude the last slice like an earlier version
% -Reads in headers first, sorts them, then reads data
% -Issues warning if SUV cannot be calculated (eg, no weight)
% -Now used dicomread_dcm2mat, which reduces the number of times that dicom
%  header information has to be read. 
% 
% Last updated: 2/2017 by Tyler Bradshaw
% Date of creation: 3/2016 by Tyler Bradshaw
% Version: 1.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parse Input %%

if length(varargin) == 1
    directory = varargin{1};   
else
    directory = uigetdir('C:\','Select directory with DICOM files');
end
%% READ IN DATA 

%loop through each file in directory (ie, slice)
list = dir(directory); %first 2 will be . and ..

zpos = zeros(1,size(list,1)-2); %axial position, for sorting
rescaleSlopes =  zeros(1,size(list,1)-2); %slopes for scaling image values
rescaleIntercepts =  zeros(1,size(list,1)-2); %intercepts for converting image values

% Get full dicominfo on a single slice to get certain info
dicomInform = dicominfo([directory '\' list(3).name]);

%Get modality and dimensions
modality = dicomInform.Modality;
Nx = double(dicomInform.Rows);
Ny = double(dicomInform.Columns);
Nz = length(list) - 2;

% preallocate
images = zeros(Nx,Ny,Nz);
info = cell(1,size(list,1)-2); %basic header info

%% loop through files, read image and zpos

% Tyler's fast method for reading images AND relevant header info at the
% same time

if strcmpi(modality, 'MR') 
    for j = 1:Nz
%         [images(:,:,j), ~, ~, ~, info{j}] = dicomread_dcm2mat([directory '\' list(j+2).name]);
        images(:,:,j) = dicomread([directory '\' list(j+2).name]);
        info{j} = dicominfo([directory '\' list(j+2).name]);
        zpos(j) = info{j}.ImagePositionPatient(3);  
    end
else

    for j = 1:Nz  
        [images(:,:,j), ~, ~, ~, info{j}] = dicomread_dcm2mat([directory '\' list(j+2).name]);
        zpos(j) = info{j}.ImagePositionPatient(3);
        rescaleSlopes(j) = info{j}.RescaleSlope;
        rescaleIntercepts(j) = info{j}.RescaleIntercept;
    end
end

[~,ordering] = sort(zpos);

%% calculate SUV, important info
    
switch modality
    case 'CT'
        %Get # z-slices and slice thickness
        CT = zeros(Nx,Ny,Nz);
        %slice thickness not in header, take different in location
        dz = abs(zpos(ordering(1)) - ...
            zpos(ordering(2))); %in mm

    case 'PT'
        %Get # z-slices and slice thickness
        if double(dicomInform.NumberOfSlices) ~= (length(list) - 2)
            warning('Check that folder only contains dicom files for this series');
        end
        dz = abs(zpos(ordering(1)) - ...
            zpos(ordering(2))); %in mm
        Bq = zeros(Nx,Ny,Nz);
        SUV = zeros(Nx,Ny,Nz);

        scan_datetime = ...
            [dicomInform.SeriesDate dicomInform.SeriesTime];             
        %if someone didn't enter the tracer info, then the inj_dose
        %will be blank. In that case, the script will write to
        %Bq/cc
        try
            %SUV factor calculation
            weight = dicomInform.PatientWeight;        %in kg
            half_life = dicomInform.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife;            %in s          
            inj_datetime = dicomInform.RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartDatetime;
            inj_dose = dicomInform.RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose;            %in Bq
            inj_datetimeVect = ...
                datevec(inj_datetime, 'yyyymmddHHMMSS.FFF');
            scan_datetimeVect = ...
                datevec(scan_datetime, 'yyyymmddHHMMSS');
            %Delay time between injection and start of scan (in s)
            t=etime(scan_datetimeVect, inj_datetimeVect);
            scan_dose = inj_dose*exp(-(log(2)*t)/half_life);%in  Bq                                            
            %Calculate SUV factor
            SUV_factor = 1/((scan_dose/weight)*0.001);
                        
             %if the SUV cannot be calucated (eg, there is no weight)
            if SUV_factor == 0
                warning('SUV cannot be calculated for this series. Units are Bq/cc');
                SUV = Bq;
                SUV_factor = 1;
            end

        catch err
            warning(['Tracer activity not found in DICOM' ...
                ' header. PET image will be in units of Bq/cc']);
            inj_dose = [];
            SUV_factor = 1;
        end
    case 'MR'
        MR = zeros(Nx,Ny,Nz);
        dz = abs(zpos(ordering(1)) - ...
            zpos(ordering(2)));
end


%% rescale images, convert to SUV
for k = 1:Nz    
    %factor to convert to correct units
    m = rescaleSlopes(ordering(k));
    b = rescaleIntercepts(ordering(k));
    
    %read data and apply scaling factors + shift
    switch modality
        case 'CT'                               
            CT(:,:,k) = double(images(:,:,ordering(k))).*m+b;
        case 'PT'         
            Bq(:,:,k) = double(images(:,:,ordering(k))).*m+b;
            SUV(:,:,k) = SUV_factor * Bq(:,:,k);
        case 'MR'
            MR(:,:,k) = double(images(:,:,ordering(k)));
    end
end

%% more info

%Number of slices, pixel size
dx = dicomInform.PixelSpacing(1); %in mm
dy = dicomInform.PixelSpacing(2); %in mm
x_start = dicomInform.ImagePositionPatient(1);
y_start = dicomInform.ImagePositionPatient(2);
z_start = zpos(ordering(1));




%% OUTPUT
Geometry = [];
header = dicomInform;
switch modality
    case 'CT'
        Geometry.data = CT;
        Geometry.voxel_size = [dx, dy, dz];
        Geometry.start = [x_start, y_start, z_start];
    case 'PT'
        if ~isempty(inj_dose) 
            Geometry.data = SUV;
            Geometry.SUV_factor = SUV_factor;
        else %in case someone didn't enter the dose
            Geometry.data = Bq;
            Geometry.SUV_factor = 1;
        end
        Geometry.voxel_size = [dx, dy, dz];
        Geometry.start = [x_start, y_start, z_start];
    case 'MR'
        Geometry.data = MR;
        Geometry.voxel_size = [dx, dy, dz];
        Geometry.start = [x_start, y_start, z_start];
          
end
