%%
clear all; close all; clc

home_fld = '/Users/ashish/Google_Drive/PhD/Mesh_Fun/Heart_Model/Biological/reSQUEEZ_Clinical/';

patient = 'cvc2007301829';

% Folder where image is stored
img_path = ['/Users/ashish/Documents/PhD/CT_Data/CRT/',patient,'/Input/img-dcm/img_00/'];

cd(img_path)
dcm_files = dir('*.dcm');

% Name of mid axial slice - must be a dicom file %
slc = 80;
if isnan(slc)
    slc = round(size(dcm_files,1)/2);
end

info = dicominfo(dcm_files(slc).name);
img = dicomread(dcm_files(slc).name).*info.RescaleSlope+info.RescaleIntercept;

%%%%%%%%%%%% Draw boundary around Myocardium + LV blood pool %%%%%%%%%%%%
figure('pos',[10 10 1200 1200]);
imagesc(img); axis equal; colormap gray; caxis([-100 400])
title('Draw Boundary Around Myocardium & LV Blood pool','FontSize',30)

h = drawpolyline;
bw = poly2mask(h.Position(:,1),h.Position(:,2),size(img,1),size(img,2));
mask = bw.*double(img);
close;

thresh = multithresh(mask,1);
thresh = round(thresh,-1);

figure('pos',[10 10 1200 1200]);
seg = zeros(size(img));
seg(mask>=thresh) = 1;
imagesc(seg); axis equal; colormap gray
title(['Threshold is: ',num2str(thresh)],'FontSize',30)

cd(home_fld)

%% Threshold from reformatted AUH data

home_fld = '/Users/ashish/Google_Drive/PhD/Mesh_Fun/Heart_Model/Biological/reSQUEEZ_Clinical/';

% Folder where image is stored
img_path = ['/Users/ashish/Documents/PhD/CT_Data/AUH/3/img-nii/3/163pre/163Pre/Unnamed - 0/lax_2ch_GLS_8747/'];

cd(img_path)
dcm_files = dir('*.dcm');

% Name of mid axial slice - must be a dicom file %
slc = 1;
if isnan(slc)
    slc = round(size(dcm_files,1)/2);
end

info = dicominfo(dcm_files(slc).name);
img = dicomread(dcm_files(slc).name).*info.RescaleSlope+info.RescaleIntercept;

%%%%%%%%%%%% Draw boundary around Myocardium + LV blood pool %%%%%%%%%%%%
figure('pos',[10 10 1200 1200]);
imagesc(img); axis equal; colormap gray; caxis([-100 400])
title('Draw Boundary Around Myocardium & LV Blood pool','FontSize',30)

h = drawpolyline;
bw = poly2mask(h.Position(:,1),h.Position(:,2),size(img,1),size(img,2));
mask = bw.*double(img);
close

thresh = multithresh(mask,1);
thresh = round(thresh,-1);

figure('pos',[10 10 1200 1200]);
seg = zeros(size(img));
seg(mask>=thresh) = 1;
imagesc(seg); axis equal; colormap gray
title(['Threshold is: ',num2str(thresh)],'FontSize',30)

cd(home_fld)


%% nii file

clear all; close all;

addpath('/Users/ashish/Google_Drive/PhD/nii_reading') 

home_fld = '/Users/ashish/Documents/PhD/CT_Data/Fractal_paper_correlation/';

patient = 'Patient1';

% Folder where image is stored
img_path = [home_fld,patient,'/Repro/Trial2/img-nii/img_00.nii.gz'];

% Reading nii file
data = load_nii(img_path);
res = data.hdr.dime.pixdim(2:4);
I = data.img;

I = permute(I,[2 1 3]);
I = flip(I,1);
I = flip(I,2);

% Slice
slc = round(size(I,3)/2);
% slc = 100;
img = I(:,:,slc);

%%%%%%%%%%%% Draw boundary around Myocardium + LV blood pool %%%%%%%%%%%%
figure('pos',[10 10 1200 1200]);
imagesc(img); axis equal; colormap gray; caxis([-200 400])
title('Draw Boundary Around Myocardium & LV Blood pool','FontSize',30)

h = drawpolyline;
bw = poly2mask(h.Position(:,1),h.Position(:,2),size(img,1),size(img,2));
mask = bw.*double(img);
close

thresh = multithresh(mask,1);
thresh = round(thresh,-1);

figure('pos',[10 10 1200 1200]);
seg = zeros(size(img));
seg(mask>=thresh) = 1;
imagesc(seg); axis equal; colormap gray
title(['Threshold is: ',num2str(thresh)],'FontSize',30)

% cd(home_fld)
