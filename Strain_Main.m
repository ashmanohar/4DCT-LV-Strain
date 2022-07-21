% Copyright (C) 2021-2022 Ashish Manohar (asmanoha@ucsd.edu)
% also see https://github.com/ashmanohar
%
%     This file is the main script to perform 4DCT-derived regional LV function analysis.
%
%     The source code is provided under the terms of the GNU General Public License as published by
%     the Free Software Foundation version 2 of the License.
% 
%     This 4DCT sytrain analysis software is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with 4DCT strain package; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA


%%% Need Statistics and Machine Learning Toolbox + Optimization Toolbox

%%% Segmentation files should be named as follows 'seg_xx' starting the count from '00'.
%   They should be located under the "patient folder" in 'seg-nii' %%%
%%% Need one img-nii file in 'img-nii' folder (under "patient folder") labeled as 'img_xx', where 'xx' corresponds to info.template %%%

%%% Need the "iso2mesh" toolbox from
%%% "http://iso2mesh.sourceforge.net/cgi-bin/index.cgi"
addpath(genpath('/Users/ashish/Google_Drive/PhD/Mesh_Fun/iso2mesh_mac'))
                              
%%% Need the nii toolbox from
% "https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image"
addpath('/Users/ashish/Google_Drive/PhD/nii_reading')

%%% Need the Coherent Point Drift poi t set registration algorithm from
% "https://sites.google.com/site/myronenko/research/cpd"
addpath(genpath('/Users/ashish/Google_Drive/PhD/Mesh_Fun/squeez-code/tools/cpd2'));        %CPD path

% Data set name
info.patient = 'cvc1910181441';

% Paths
info.home_path = '/Users/ashish/Google_Drive/PhD/4DCT_Strain/';   %Home folder containing function scripts
addpath(genpath(info.home_path))
info.save_path = ['/Users/ashish/Google_Drive/PhD/Test_Software/',...
    info.patient,'/'];     %Folder path to save results
mkdir(info.save_path)
info.seg_path = ['/Users/ashish/Documents/PhD/CT_Data/ForTrainingMovie/',...
    info.patient,'/Input/seg-nii/'];       %Folder path containing the segmentations to read
info.img_path = ['/Users/ashish/Documents/PhD/CT_Data/ForTrainingMovie/',...
    info.patient,'/Input/img-nii/'];       %Folder path containing the nii images to read for rotations
info.dcm_path = ['/Users/ashish/Documents/PhD/CT_Data/ForTrainingMovie/',...
    info.patient,'/Input/img-dcm/'];       %Folder path containing the dcm images to read for percent rr

% Image preparation parmeters
info.iso_res = 0.5;                 %Standardized starting resolution in mm
info.desired_res = 2;               %Desired operating resolution in mm (multiples of 0.5)
info.averaging_threshold = 0.5;     %Threshold for voxels post averaging
info.fill_paps = 0;                 %Flag for filling in pap muscles using 3D convex hull: Need MATLAB 2017b or above

% Segmentation labels
info.lv_label = 1;
info.la_label = 2;
info.lvot_label = 4;

% Finding percent RR and time between frames - Works only for GE scans. For
% other manufacturers, have to change the corresponding dicom tags to read
% gantry angle and cardiac phase percentage.
%%% Need a minimum of info.percent_rr to make strain vs %R-R curves
info = Timing_Info(info);
disp(['Percent R-R: ',num2str(info.percent_rr)])
disp(['Gantry difference: ',num2str(round(diff(info.gantry_angle)))])
disp(['Number of timeframes: ',num2str(length(info.percent_rr))])

% Time frame list
info.timeframes = 1:length(info.percent_rr);             %Desired time frames for analysis
info.template = 1;                  %Used as template mesh for registration
info.reference = 1;                 %Used as reference phase for strain calculations

%Flag for computing endocardial circ. and longitudinal strains
info.endo_strains = 1;

str = input(['Time frames ',num2str(info.timeframes(1)),' to ',num2str(info.timeframes(end)),...
    ' and template/reference frame of ',num2str(info.template),' ok? y/n: '],'s');
if str ~= 'y'
    error('Please check time frames');
end

%Flag for temporally smoothing CPD vertices using Fourier decomposition.
%Use only when there is one complete period (1 full cardiac cycle)
info.smooth_verts = 1;
if info.smooth_verts == 1 && length(info.timeframes) <= 5
    error('Please check temporal smoothing flag for periodicity of entered time frames')
elseif info.smooth_verts == 0 && length(info.timeframes) >= 5
    error('Temporal smoothing switched off')
end

% Axes limits for plotting
%When plotting use 'xlim([info.xlim]); ylim([info.ylim]); zlim([info.zlim])' for non-rotated LV
%When plotting use 'xlim([info.rot_xlim]); ylim([info.rot_ylim]); %zlim([info.rot_zlim])' for rotated LV

clear str
disp('xxxxxxxxx - Analysis Parameters Saved - xxxxxxxxx')


%% Image Processing and Mesh Extraction 

[Mesh,info] = Mesh_Extraction(info);

disp('xxxxxxxxx - Meshes Extracted - xxxxxxxxx')


%% Registration - SQUEEZ

% CPD Parameters
opts.corresp = 1;
opts.normalize = 1;
opts.max_it = 1500;
opts.tol = 1e-5;
opts.viz = 0;
opts.method = 'nonrigid_lowrank';
opts.fgt = 2;
opts.eigfgt = 0;
opts.numeig = 100;
opts.outliers = 0.05;
opts.beta = 2;
opts.lambda = 3;

Mesh = Registration(Mesh,info,opts);

disp('xxxxxxxxx - Registration done - xxxxxxxxx')

clear opts


%% Mesh Rotation

[Mesh, info] = Rotation(Mesh,info);

disp('xxxxxxxxx - Meshes rotated - xxxxxxxxx')


%% Calculating Regional Shortening

Mesh = RSCT(Mesh,info);

disp('xxxxxxxxx - Regional shortening calculated - xxxxxxxxx')


%% Calculating Regional Circumferential and Regional Longitudinal Strain

if info.endo_strains
    Mesh = Ecc_Ell(Mesh,info);
    disp('xxxxxxxxx - Regional endocardial strains calculated - xxxxxxxxx')
else
    disp('xxxxxxxxx - Regional endocardial strains NOT calculated - xxxxxxxxx')
end


%% Polar Sampling of regional shortening and endocardial strains
% Creating raw data set of "high resolution" sampling of strain values as a function of theta and z

info.rawdata_slicethickness = 5/info.desired_res;           %Enter slice thickness for raw data sampling in mm. Has to be >2*info.desired_res
info.apical_basal_threshold = [0.05 0.05];                  %Apical and basal percentage tolerance in that order

if info.rawdata_slicethickness <= 2
    error('Slice thickness too small')
else
    [Mesh, info] = Data_Sampling(Mesh,info);
    disp('xxxxxxxxx - Polar data sampled - xxxxxxxxx')
end


%% AHA Plotting
% Hard coded 16 AHA segments

info.RSct_limits = [-0.5 0.1];

%If the "Timing_Info" did not work, manually input the cardiac phase
%percentages here
% info.percent_rr = [0:5:95];

Mesh = AHA(Mesh,info);

disp('xxxxxxxxx - AHA plots generated - xxxxxxxxx')


%% Bullseye Plotting

info.polar_res = [36 10];                       %Enter desired number of points in bullseye plots in the format number [azimuthal radial]
info.polar_NoOfCols = 5;                        %Number of columns in bullseye plot subplot
info.RSct_limits = [-0.3 0.1];

Mesh = Bullseye_Plots(Mesh,info);

disp('xxxxxxxxx - Bullseye plots generated - xxxxxxxxx')


%% High-resolution strain maps - 90 segments; 5 z slices and 18 in theta

%desired resolution in [azimuthal direction, slice direction]
info.aha_highres = [18 5]; % do not recommend going higher than this for 2mm resolution meshes

info.RSct_limits = [-0.5 0.1];

[Mesh, info] = AHA_Highres(Mesh,info);

disp('xxxxxxxxx - Hi-Res maps generated - xxxxxxxxx')


%% Saving variables

save([info.save_path,info.patient,'.mat'],'Mesh','info')

disp('xxxxxxxxx - Done saving variables - xxxxxxxxx')

