clear
orig_folder = pwd;
addpath(genpath([pwd]))
cd('E:\resarch_data\fMRI\BrainVortexToolbox-main')
load("parcellation_template.mat")
load("parcellation_template22_RightBrain_subject1-100.mat")
hemisphere = 1; % 1 for left hemisphere, 2 for right hemisphere, keep this as 1 for demonstration purpose


% a range of parameters avaiable for different dataset, but for demonstration
% purpose, only use the parameters provided


timelimit = 216; 

sbjlist = 1:20;
% func = 'DMT';
func = 'PCB';

%% MAIN SUBJECT FOR LOOP

% For each time step, there will be: x, y, color, size, opacity
% initialise a structure to store all the data
% world = [];
% for i = 1:timelimit
%     world_slices = struct('x', [], 'y', [], 'spin', [], 'size', [], 'opacity', []);
%     world = [world, world_slices];
% end



all_spirals_neg = {};
all_spirals_pos = {};
all_spirals_overall = {};

for i = 1:size(sbjlist,2)
% for i = 1:1
    fprintf("Processing for sub" + i + "\n")
    tic

    if func == 'DMT'
        filename = ['E:\resarch_data\fMRI\fMRI_DMT\BrainToolBoxed\DMT_post1\Spiral Detected\', 'Spiral_detected_surfilt_motor_task_LEFT_sub',sprintf('%d',i),'.mat']; 
    elseif func == 'PCB'
        filename = ['E:\resarch_data\fMRI\fMRI_DMT\BrainToolBoxed\PCB_post1\Spiral Detected\', 'Spiral_detected_surfilt_motor_task_LEFT_sub',sprintf('%d',i),'.mat']; 
    end
    load(filename, 'spiral_filt_nega_real_95perc_extend', 'spiral_filt_pos_real_95perc_extend')
    toc



    spiral_filt_nega = spiral_filt_nega_real_95perc_extend; % Spiral singularity locations
    spiral_filt_pos = spiral_filt_pos_real_95perc_extend;
    
%    existing_period = zeros(size(spiral_filt_nega));


    % anticlockwise spirals
    braingrid_pos = zeros(175, 251);
    braingrid_neg = zeros(175, 251);
    braingrid_overall = zeros(175, 251);

    for ipatt = 1:size(spiral_filt_nega,1)
        for t = 1:size(spiral_filt_nega,2)
            if ~isempty(spiral_filt_nega{ipatt,t})
                braingrid_neg = braingrid_neg + spiral_filt_nega{ipatt,t};
                braingrid_overall = braingrid_overall + spiral_filt_nega{ipatt,t};
            end
        end
    end

    % clockwise spirals
    for ipatt = 1:size(spiral_filt_pos,1)
        for t = 1:size(spiral_filt_pos,2)
            if ~isempty(spiral_filt_pos{ipatt,t})
                braingrid_pos = braingrid_pos + spiral_filt_pos{ipatt,t};
                braingrid_overall = braingrid_overall + spiral_filt_pos{ipatt,t};
            end
        end
    end

    all_spirals_neg{i} = braingrid_neg;
    all_spirals_pos{i} = braingrid_pos;
    all_spirals_overall{i} = braingrid_overall;
end



if func == 'DMT'
    fn = ['E:\resarch_data\fMRI\fMRI_DMT\BrainToolBoxed\DMT_post1\forplotting\','alltime_spiral_grid_DMT.mat']
    save(fn, 'all_spirals_neg', 'all_spirals_pos', 'all_spirals_overall')
elseif func == 'PCB'
    fn = ['E:\resarch_data\fMRI\fMRI_DMT\BrainToolBoxed\PCB_post1\forplotting\','alltime_spiral_grid_PCB.mat']
    save(fn, 'all_spirals_neg', 'all_spirals_pos', 'all_spirals_overall')
end

cd(orig_folder)