% dataDir = 'E:\resarch_data\fMRI\fMRI_Motor_Data\Raw_Data\102109\MNINonLinear\Results\tfMRI_MOTOR_LR\tfMRI_MOTOR_LR_Atlas.dtseries.nii';
addpath(genpath([pwd]))
dataDir = 'E:\resarch_data\fMRI\ciftiy_test\out\sub02\MNINonLinear\Results\LSD_Rest1_clean_test2\LSD_Rest1_clean_test2_Atlas_s0.dtseries.nii';
posFile = 'E:\resarch_data\fMRI\fMRI_Motor_Data\Data_Pos\102109\MNINonLinear\fsaverage_LR32k\102109.L.flat.32k_fs_LR.surf.gii';

load("parcellation_template.mat")


flagSur = 0;

params.sigmScale =  [29.35 14.93];            % bandpass 5 bandwidth ranges
params.downSRate = 2 ;                        % downsample the re-interpolation                  
params.xCord = -250:params.downSRate:250 ;    % coordinate re-interpolation, left hemisphere
params.yCord = -150:params.downSRate:200 ;    
params.fsTem = 1/0.72 ;                       % temporal sampling rate


flagRest = 0;
flagSmooth = 2;
hemisphere = 1;
flagTask = 4;



disp('start importing and pre-processing ...')

downSRate = params.downSRate ;
xCord =params.xCord ;
yCord = params.yCord ;

%% Import cortex data
tic
% addpath(genpath([pwd,'/ToolOthers/fMRI/cifti-matlab-master']))
% load the surface data
data = ft_read_cifti(dataDir) ;
% load the position data
posLC = gifti(posFile) ;

fsTem = 1/0.72 ;
% this function extract the cortex part of data from the HCP data
[sigValid,posValid,nanChans] = preproc_fRMI(data,posLC,fsTem) ;
toc
% clearvars data posLC posRC

% checkRegion(posValid)

%% new interpolation based on the whole left cortex
x = double(posValid.vertices(:,1));
y = double(posValid.vertices(:,2));
% k = convhull(x,y);
k = alphaShape(x,y,4) ;
[a, b] = k.boundaryFacets();

bw = poly2mask(b(:,1)-min(xCord)+1,b(:,2)-min(yCord)+1,max(yCord)-min(yCord)+1,...
    max(xCord)-min(xCord)+1);

flagPlotSpe = 0 ;
flagCheckIntp = 0 ;

mask = double(bw(1:downSRate:end,1:downSRate:end)) ;
mask(mask==0) = nan ;
data_reInterp = spaceFreq_fMRI(posValid,sigValid,xCord,yCord,...
    flagPlotSpe, flagCheckIntp,mask) ;
dataOut = data_reInterp;



%% temporal bandpass filtering

dataIn = dataOut;
fLow = 0.01 ; % lower limit of the bandpass filter
fHigh = 0.1 ; % upper limit of the bandpass filter
dataOut_reshape = reshape(dataIn,size(dataIn,1)*size(dataIn,2),[]) ;
[bandpasSig,~, ~, ~] = bandpa_fMRI(dataOut_reshape,fsTem,fLow,fHigh) ;
%     bandpasSig = zscore(bandpasSig,[],2) ; % ***************
bandpass_reInterp = reshape(bandpasSig,size(dataIn,1),size(dataIn,2),[]) ;
dataOut = bandpass_reInterp;


%% spatial bandpass filtering - Gaussian

dataIn = dataOut(1:175,:,:);
sigmScale = params.sigmScale/downSRate ;
sigLPass = [] ;
sigBPass = [] ;
% spatial bandpass filter
numScale = length(sigmScale)-1 ;
for iTime = 1:size(dataIn,3)
    filtSigma = sigmScale(1);   % 0.6
    filtWidth = ceil(3*filtSigma);
    imageFilter=fspecial('gaussian',filtWidth,filtSigma);
    sigLPass(1,:,:,iTime) = nanconv(dataIn(:,:,iTime),imageFilter,'edge', 'nanout');
    for iScale = 1:numScale
        filtSigma = sigmScale(iScale+1);
        filtWidth = ceil(3*filtSigma);
        imageFilter=fspecial('gaussian',filtWidth,filtSigma);
        sigLPass(iScale+1,:,:,iTime) = nanconv(dataIn(:,:,iTime),imageFilter,'edge', 'nanout');
        sigBPass(iScale,:,:,iTime) = sigLPass(iScale+1,:,:,iTime) - sigLPass(iScale,:,:,iTime) ;
    end
end

dataOut = permute(sigBPass,[2,3,4,1]);
dataOut = dataOut(1:175,:,:).*parcellation_template(1:175,:)./parcellation_template(1:175,:);




timelimit = 217; % For LSD Data
% timelimit = 283; % For Motor Data

DataIn_smooth = dataOut;


% Code for calculating phase vector field
phaseSig = DataIn_smooth;
vPhaseX = zeros(size(phaseSig)) ;
vPhaseY = zeros(size(phaseSig)) ;
% cd(main_folder)
for iTime = 1:size(phaseSig,3)
    for iX = 1:size(phaseSig,1)
        vPhaseX(iX,2:end-1,iTime) = (anglesubtract(phaseSig(iX,3:end,iTime),phaseSig(iX,1:end-2,iTime)))/2 ;
    end
    for iY = 1:size(phaseSig,2)
        vPhaseY(2:end-1,iY,iTime) = (anglesubtract(phaseSig(3:end,iY,iTime),phaseSig(1:end-2,iY,iTime)))/2 ;
    end
end
Vx_flowmap_norm_phase = -vPhaseX./sqrt(vPhaseX.^2+vPhaseY.^2) ; % standardize vector length as 1
Vy_flowmap_norm_phase = -vPhaseY./sqrt(vPhaseX.^2+vPhaseY.^2) ; % standardize vector length as 1

%% To Visulaise


time = 20;


% Phase Map Plotting
figure(2)

imagesc(DataIn_smooth(:,:,time))
colorbar()
set(gca,'ydir','normal')
title('Phase Map')
axis equal
set(gcf,  'unit', 'centimeters', 'position', [15 12 26 20]); 


% Phase Vector field plotting
figure(3)

hold on
quiver(Vx_flowmap_norm_phase(:,:,time),Vy_flowmap_norm_phase(:,:,time)) 
for parcellation_ID = 1:22 % Numbers correspond to certain areas of brain
    parcellation_template_1par = parcellation_template;
    parcellation_template_1par(isnan(parcellation_template_1par)) = 0;

    parcellation_template_1par(parcellation_template_1par~=parcellation_ID) = 0;

    B = bwboundaries(parcellation_template_1par(:,:,1),'noholes');
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2),boundary(:,1),'-','linewidth',2,'color',[0,0,0])

    end
end
hold off
title('Phase Vector field plotting')
axis equal
set(gcf,  'unit', 'centimeters', 'position', [15 12 26 20]); 




