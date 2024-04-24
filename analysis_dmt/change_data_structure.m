clear

% read data
filepaths = get_spiral_detected_path("LSD", "LSD");

for i = 1:1
    spiraldetecetd = load(filepaths{i});
end