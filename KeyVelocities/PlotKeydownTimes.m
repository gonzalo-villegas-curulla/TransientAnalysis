files = dir('*.mat');

for idx = 1 : length(files)
    
    load(files(idx).name);
    disp(length(MX));
end