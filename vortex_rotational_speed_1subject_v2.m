function [phase_vector_angle_speed_pos_nega] = vortex_rotational_speed_1subject_v2(subject,flagRest)


%%  Vortex rotational speed
disp(['initiating process...'])    
disp(['initating processing of vortex rotational speed for single subject...'])
disp(['starting subject ',num2str(subject)])    

%% load vortex centre data
disp(['loading preprocessed fMRI data files, vortex centre coordinates and full-sized vortex distributions...'])    

cd '/import/headnode2/yixu4976/Downloads/fMRI/PBS_access/Vortex Detected/Combatibility_ratio_based_surrogate_filtered/Bandpass_29.35_14.93_Rest1199_30DegreePhaseDifferenceThreshold_95percfilt_SignificantFrameExtend_WithinVortex'
filename = ['vortex_bandpass_29.35_14.93_Rest1199_30DegPhaDif_compatibility_95perc_SignificantFrameExtend_sub',num2str(subject),'.mat'];
% cd '/import/headnode2/yixu4976/Downloads/fMRI/PBS_access/Vortex Detected/Bandpass_Resting1200_LEFT_29.35_14.93_1sThreshold_FDRadj_FullsizePLusCentre'
% filename = ['vortex_LEFT_bandpass_29.35_14.93_1sThreshold_FDRadj_FullsizePlusCentre_subject',num2str(subject),'.mat'];
load(filename)      
vortex_filt_nega = vortex_filt_nega_real_FDRadj_extend;
vortex_filt_pos = vortex_filt_pos_real_FDRadj_extend;    
vortex_filt_nega_centre = vortex_filt_nega_real_centreONLY_FDRadj_extend;
vortex_filt_pos_centre = vortex_filt_pos_real_centreONLY_FDRadj_extend;
%     vortex_filt_nega = vortex_filt_nega_maxsize_99perc_FDRadj;
%     vortex_filt_pos = vortex_filt_pos_maxsize_99perc_FDRadj;    
%     vortex_filt_nega_centre = vortex_filt_nega_maxsize_99perc_FDRadj_centreONLY;
%     vortex_filt_pos_centre = vortex_filt_pos_maxsize_99perc_FDRadj_centreONLY;

clearvars -except vortex_filt_nega vortex_filt_pos vortex_filt_nega_centre vortex_filt_pos_centre subject
% load fMRI amplitude data
%     filename = ['Preprocessed_resting_bandpass_ScaleSigma_29.35_14.93_sub',num2str(subject),'.mat'];
cd '/import/headnode2/yixu4976/Downloads/fMRI/Data_Raw/toolbox/preprocess/Preprocessed_resting_LeftBrain_bandpass_ScaleSigma_29.35_14.93_v2'
filename = ['Preprocessed_resting_LEFT_bandpass_ScaleSigma_29.35_14.93_sub',num2str(subject),'_v2.mat'];    
load(filename)     
    


%%  Find phase map of fMRI signal
cd '/import/headnode2/yixu4976/Downloads/fMRI/Manuscript/Code For Github_v2'
phase_map = zeros(176,251,1200);
for irow = 1:176
    for icol = 1:251 
        temp1_sig = reshape(sigBPass(:,irow,icol,:),[1200,1]);
        if nansum(temp1_sig(:)) ~= 0
            phase_map(irow,icol,:) = angle(hilbert(temp1_sig));
        end 
    end
end


% Filter phase map within vortex
% phase_map_filt_nega = [];
% phase_map_filt_pos = [];
phase_vector_angle_speed_nega = [];
phase_vector_angle_speed_pos = [];

% negative (clockwise) vortex only

for ipatt = 1:size(vortex_filt_nega,1)
    phase_vector_angle_1patt = [];
%     phase_map_filt_1t = zeros(176,251);
    for t = 1:size(vortex_filt_nega,2)
        temp1 = full(vortex_filt_nega{ipatt,t});
        vortex_nega_1t1patt = zeros(176,251);
        if nansum(temp1(:)) ~= 0
%            vortex_nega_1t1patt(1:175,1:251) = temp1;
           vortex_nega_1t1patt(1:176,1:251) = temp1;
           vortex_filt_nega_centre_1t1patt = vortex_filt_nega_centre{ipatt,t};
           phase_map_filt_1t1patt = abs(vortex_nega_1t1patt).*phase_map(:,:,t);
           phase_map_filt_1t1patt_pifilt = phase_map_filt_1t1patt;
           phase_map_filt_1t1patt_pifilt_abs = abs(phase_map_filt_1t1patt);
           phase_map_filt_1t1patt_pifilt(phase_map_filt_1t1patt_pifilt_abs<5/8*pi) = 0;           
           % vectorize the angle between ~Pi phase points and vortex centers
           idx = find(phase_map_filt_1t1patt_pifilt);
           if isempty(idx) == 1
               phase_vector_angle_1patt = [phase_vector_angle_1patt; nan];
           end
           [row col] = ind2sub([176 251],idx);
           phase_vector_u_nega = nanmean(col-vortex_filt_nega_centre_1t1patt(1));
           phase_vector_v_nega = nanmean(row-vortex_filt_nega_centre_1t1patt(2)); 
           phase_vector_angle_1t1patt = angle(phase_vector_u_nega + i*phase_vector_v_nega);
           phase_vector_angle_1patt = [phase_vector_angle_1patt ; phase_vector_angle_1t1patt];
        end  
    end
    if nansum(phase_vector_angle_1patt(:)./phase_vector_angle_1patt(:)) > 1
       phase_vector_angle_1patt_diff = anglesubtract(phase_vector_angle_1patt(2:end),phase_vector_angle_1patt(1:end-1));
       phase_vector_angle_speed_nega = [phase_vector_angle_speed_nega; phase_vector_angle_1patt_diff];
    end  
end
for i2 = 1:size(phase_vector_angle_speed_nega,1)
    if phase_vector_angle_speed_nega(i2) > pi
        phase_vector_angle_speed_nega(i2) = phase_vector_angle_speed_nega(i2) - 2*pi;
    elseif phase_vector_angle_speed_nega(i2) < -pi
        phase_vector_angle_speed_nega(i2) = phase_vector_angle_speed_nega(i2) + 2*pi;
    end
end
% Positive (anti-clockwise) vortex only

for ipatt = 1:size(vortex_filt_pos,1)
    phase_vector_angle_1patt = [];
%     phase_map_filt_1t = zeros(176,251);
    for t = 1:size(vortex_filt_pos,2)
        temp1 = full(vortex_filt_pos{ipatt,t});
        vortex_pos_1t1patt = zeros(176,251);        
        if nansum(temp1(:)) ~= 0
%            vortex_pos_1t1patt(1:175,1:251) = temp1;  
           vortex_pos_1t1patt(1:176,1:251) = temp1;
           vortex_filt_pos_centre_1t1patt = vortex_filt_pos_centre{ipatt,t};
           phase_map_filt_1t1patt = abs(vortex_pos_1t1patt).*phase_map(:,:,t);
           phase_map_filt_1t1patt_pifilt = phase_map_filt_1t1patt;
           phase_map_filt_1t1patt_pifilt_abs = abs(phase_map_filt_1t1patt);
           phase_map_filt_1t1patt_pifilt(phase_map_filt_1t1patt_pifilt_abs<5/8*pi) = 0;           
           % vectorize the angle between ~Pi phase points and vortex centers
           idx = find(phase_map_filt_1t1patt_pifilt);
           if isempty(idx) == 1
               phase_vector_angle_1patt = [phase_vector_angle_1patt; nan];
           end
           [row col] = ind2sub([176 251],idx);
           phase_vector_u_pos = nanmean(col-vortex_filt_pos_centre_1t1patt(1));
           phase_vector_v_pos = nanmean(row-vortex_filt_pos_centre_1t1patt(2)); 
           phase_vector_angle_1t1patt = angle(phase_vector_u_pos + i*phase_vector_v_pos);
           phase_vector_angle_1patt = [phase_vector_angle_1patt ; phase_vector_angle_1t1patt];
        end  
    end
    if nansum(phase_vector_angle_1patt(:)./phase_vector_angle_1patt(:)) > 1
       phase_vector_angle_1patt_diff = anglesubtract(phase_vector_angle_1patt(2:end),phase_vector_angle_1patt(1:end-1));
       phase_vector_angle_speed_pos = [phase_vector_angle_speed_pos; phase_vector_angle_1patt_diff];
    end
end

for i2 = 1:size(phase_vector_angle_speed_pos,1)
    if phase_vector_angle_speed_pos(i2) > pi
        phase_vector_angle_speed_pos(i2) = phase_vector_angle_speed_pos(i2) - 2*pi;
    elseif phase_vector_angle_speed_pos(i2) < -pi
        phase_vector_angle_speed_pos(i2) = phase_vector_angle_speed_pos(i2) + 2*pi;
    end
end

phase_vector_angle_speed_pos_nega = abs([phase_vector_angle_speed_pos(:);phase_vector_angle_speed_nega(:)]);

%  phase_vector_angle_speed_pos_avg = nanmean(abs(phase_vector_angle_speed_pos(:)));
%  phase_vector_angle_speed_nega_avg = nanmean(abs(phase_vector_angle_speed_nega(:)));

 
filename2 = ['vortex_rotational_speed_rest1199_LEFT_sub',num2str(subject),'.mat'];
% save(['/import/headnode2/yixu4976/Downloads/fMRI/Manuscript/Revision/vortex_statistics/',filename2],'phase_vector_angle_speed_pos_nega','subject','-v7.3')

end


