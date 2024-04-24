function paths = get_spiral_detected_path(func, group)

    if strcmpi(func, "LSD")
        datapath = "E:\resarch_data\fMRI\fMRI_LSD\Spiral Detected";
        fnformat = "Spiral_detected_surfilt_%s_LEFT_S%02d.mat";

    end

    if strcmpi(group, "LSD")
        group = "LSD";
        id = 1:14;
    elseif strcmpi(group, "PCB")
        group = "PCB";
        id = 1:14;
    end

    paths = arrayfun(@(x) ...
        fullfile(datapath,sprintf(fnformat, group, x)), ...
        id, 'UniformOutput', false);

    
end