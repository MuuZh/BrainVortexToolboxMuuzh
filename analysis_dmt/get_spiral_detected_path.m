function paths = get_spiral_detected_path(func, group)

    keys =["LSD_LSD", "LSD_PCB", "DMT_DMT", "DMT_PCB"];

    data_paths = ["E:\resarch_data\fMRI\fMRI_LSD\Spiral Detected";
                "E:\resarch_data\fMRI\fMRI_LSD\Spiral Detected";
                "E:\resarch_data\fMRI\fMRI_DMT\BrainToolBoxed\DMT_post1\Spiral Detected";
                "E:\resarch_data\fMRI\fMRI_DMT\BrainToolBoxed\PCB_post1\Spiral Detected"
                ];

    fn_formats = ["Spiral_detected_surfilt_LSD_LEFT_S%02d.mat"; 
                "Spiral_detected_surfilt_PCB_LEFT_S%02d.mat";
                "Spiral_detected_DMTpost1_sub%02d.mat";
                "Spiral_detected_PCBpost1_sub%02d.mat"
                ];

    data_ids = {1:14; 1:14 ; 1:22; 1:20};

    valueStructs = cellfun(@(path, format, ids) struct(...
    'DataPath', path, ...
    'FileNameFormat', format, ...
    'DataID', ids), ...
    data_paths, fn_formats, data_ids, 'UniformOutput', false);

    values = containers.Map(keys, valueStructs);
    
    keyToLookup = upper(sprintf('%s_%s', func, group));

    if isKey(values, keyToLookup)
        data = values(keyToLookup);
        if ~isempty(data.DataPath)
            % fprintf('Data Path: %s\n', data.DataPath);
            % fprintf('File Name Format: %s\n', data.FileNameFormat);
            % fprintf('Data ID Range: %s\n', mat2str(data.DataID));

            paths = arrayfun(@(id) fullfile(data.DataPath, sprintf(data.FileNameFormat, id)), data.DataID, 'UniformOutput', false);
        else
            disp('No data available for this key.');
        end
    else
        disp('Key not found.');
    end

end