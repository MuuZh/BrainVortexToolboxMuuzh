function [DataOut] = Preprocessing_main(subject,hemisphere,flagSur,flagRest,flagTask,flagSmooth)
tStart = tic;
disp(['initiating process...'])    
% add all the subfunctions
addpath(genpath([pwd]))
main_folder = pwd;
task_folder_path = '/data/DMT_data/PCB_group';
working_directory_path = [main_folder, task_folder_path];
raw_data_path = [working_directory_path,'/Raw Data'];
data_pos_path = [working_directory_path,'/Data Pos'];
preprocessed_path = [working_directory_path,'/Preprocessed Data'];


raw_data_dir = dir(raw_data_path);
data_pos_dir = dir(data_pos_path);
preprocessed_dir = dir(preprocessed_path);






% set basic parameters
if flagSmooth == 0 && flagSur == 1 % no need to do surrogate raw data
    return
end
if flagSmooth == 2 && flagSur == 1 % no need to do smoothed surrogate data as it is done already when flagSmooth == 1
    DataOut = [];
    return
end
%%
params.sigmScale =  [29.35 14.93];            % bandpass 5 bandwidth ranges
params.downSRate = 2 ;                        % downsample the re-interpolation
if hemisphere == 1                    
params.xCord = -250:params.downSRate:250 ;    % coordinate re-interpolation, left hemisphere
params.yCord = -150:params.downSRate:200 ;    
elseif hemisphere == 2                
params.xCord = -270:params.downSRate:230 ;    % coordinate re-interpolation, right hemisphere
params.yCord = -180:params.downSRate:170 ;    
end
params.fsTem = 1/0.72 ;                       % temporal sampling rate

% load raw and position data files
for iSub = subject:subject
    tic
    if flagRest == 1 % resting data
        fprintf("flagRest == 1, check parameter\n")
    elseif flagRest == 0 % task data
        if flagTask == 4 % DMT data PCB group, original 100 subjects
            dataDir = [raw_data_dir(iSub+2).name]; 
            if hemisphere == 1  
                posFile = [data_pos_path,'/fsaverage_LR32k/L.flat.32k_fs_LR.surf.gii'] ; % Left brain
            elseif hemisphere == 2  
                posFile = [data_pos_path,'/fsaverage_LR32k/R.flat.32k_fs_LR.surf.gii'] ; % Right brain
            end
        end
    end
    % process the data files
    cd(main_folder)
    surMethodNum = 7 ; flagVisBpSig = 0 ;   
    DataOut = load_fMRI(dataDir,posFile,flagSur,surMethodNum,params,flagVisBpSig,flagRest,flagSmooth,hemisphere,flagTask) ;

    % real data   
    if flagSur == 0
       if flagRest == 1
            fprintf("flagRest == 1, check parameter\n")
       elseif flagRest == 0
          if flagTask == 4
            cd(preprocessed_path)
              if flagSmooth == 1
                if hemisphere == 1 % left hemisphere
                     save(['Preprocessed_temporalbandpass_data_motor_task_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     save(['Preprocessed_temporalbandpass_data_motor_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end                     
              elseif flagSmooth == 0
                  if hemisphere == 1 % left hemisphere
                     save(['Preprocessed_raw_data_motor_task_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     save(['Preprocessed_raw_data_motor_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end       
              elseif flagSmooth == 2
                  if hemisphere == 1 % left hemisphere
                     save(['Preprocessed_spatiotemporalbandpass_data_motor_task_LEFT_sub',num2str(iSub),'.mat'],'DataOut')   
                  elseif hemisphere == 2 % right hemisphere
                     save(['Preprocessed_spatiotemporalbandpass_data_motor_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut')                 
                  end                      
              end
          end
           
       end
       
    % surrogate data (null model)  
    elseif flagSur == 1
       if flagRest == 1
            fprintf("flagRest == 1, check parameter\n")
       elseif flagRest == 0
          if flagTask == 4 % working memory task
            cd(preprocessed_path)
              if flagSmooth == 1 
                  if hemisphere == 1 % left hemisphere
                     DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                     save(['Preprocessed_temporalbandpass_data_sur_motor_task_LEFT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')   
                     DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                     save(['Preprocessed_spatiotemporalbandpass_data_sur_motor_task_LEFT_sub',num2str(iSub),'.mat'],'DataOut_smooth')                 
                  elseif hemisphere == 2 % right hemisphere
                     DataOut_unsmooth = DataOut(:,:,:,2);  % unsmoothed
                     save(['Preprocessed_temporalbandpass_data_sur_motor_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_unsmooth')                 
                     DataOut_smooth = DataOut(:,:,:,1);  % smoothed
                     save(['Preprocessed_spatiotemporalbandpass_data_sur_motor_task_RIGHT_sub',num2str(iSub),'.mat'],'DataOut_smooth')                 
                  end                                              
              end
          end
           
       end       
        
        
    end
    


        disp(['finishing subject ',num2str(iSub)])
        toc
        cd(main_folder)
end
tEnd = toc(tStart)
end


