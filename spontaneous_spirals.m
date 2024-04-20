function spontaneuous_spiral_perFrame = spontaneous_spirals(subject)
%% BASELINE Spiral time distribution
restoredefaultpath
cd '/headnode2/akar5239/BrainVortexToolbox-main/'
addpath(genpath([pwd]))
main_folder = pwd;

flagSur = 0;%  0 for real data, 1 to generate surrogate data

% a range of parameters avaiable for different dataset, but for demonstration
% purpose, only use the parameters provided
No_of_Subject = 1; % number of subjects used for analysis, randomly selected from HCP database (S1200)
flagRest = 0; %  resting data => 1 , task data => 0

hemisphere = 2; % 1 for left hemisphere, 2 for right hemisphere, keep this as 1 for demonstration purpose

% flagTask:
% 1 = language task, original 100 subjects;
% 2 = language task, additional 100 subjects;
% 3 = working memory task;
% 4 = Motor Task
flagTask = 4;

% 0 = unsmoothed, raw data; 1 = temporally smoothed (bandpass filtered)
% data; 2 = spatiotemporally smoothed (bandpass filtered) data
flagSmooth = 1;


foldername = [main_folder,'/Sample Data/Motor Task/Spiral Detected'];
cd(foldername)
if hemisphere == 1
    filename = ['Spiral_detected_surfilt_motor_task_LEFT_sub',num2str(subject),'.mat'];
    load(filename)
elseif hemisphere == 2
    filename = ['Spiral_detected_surfilt_motor_task_RIGHT_sub',num2str(subject),'.mat'];
    load(filename)
end



%mean number of new spirals for each frame

spontaneuous_spiral_perFrame = zeros(1,size(spiral_filt_pos_real,2));
for time = 2:size(spiral_filt_pos_real,2)
    spontaneuous_spiral_count = 0;
    for ipatt = 1:size(spiral_filt_pos_real,1)
        if isempty(spiral_filt_pos_real{ipatt,time-1}) && ~isempty(spiral_filt_pos_real{ipatt,time})
            spontaneuous_spiral_count = spontaneuous_spiral_count + 1;
        end
    end
    for ipatt = 1:size(spiral_filt_nega_real,1)
        if isempty(spiral_filt_nega_real{ipatt,time-1}) && ~isempty(spiral_filt_nega_real{ipatt,time})
            spontaneuous_spiral_count = spontaneuous_spiral_count + 1;
        end
    end
    spontaneuous_spiral_perFrame(time) = spontaneuous_spiral_count;
end

foldername = [main_folder,'/Sample Data/Motor Task/Analysis/'];
filename = ['spontaneuous_spirals',num2str(subject),'.mat'];
save([foldername,filename],'spontaneuous_spiral_perFrame');


end
