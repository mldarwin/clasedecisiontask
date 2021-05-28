function [x] = pullData(sID)
matToGet = dir(sprintf('VNI_behavior_VNI%s_*.mat',sID)); %.mat file
x = load(matToGet(end).name); % load the .mat file for participant (subjdata)
end

