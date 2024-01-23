function [mask_path,rawohnepath,rawmitpath,trajpath,baseoutbasepath,rotatedeg] = get_config(number, path_prefix)

% set all the paths for finding and saving necessary data and reco results

if number == 1
%%%%%%%%%% Example Volunteer 1

basepath = "C://Temp/rawData/";
mask_path = strcat(basepath,"recons/mask_all.mat");
rawohnepath = strcat(basepath,"mat_files/DataFileOfUnpreparedMeasurement_Volunteer1.mat");
rawmitpath = strcat(basepath,"mat_files/DataFileOfIRpreparedMeasurement_Volunteer1.mat");
trajpath = strcat(basepath,"../traj.mat");
baseoutbasepath = strcat(path_prefix, "/recons/volunteer1_");
rotatedeg = 131; % might be needed to ensure AP is up-down in final recos

elseif number == 2
%%%%%%%%%% Example Volunteer 2

basepath = "C://Temp/rawData/";
mask_path = strcat(basepath,"recons/mask_all.mat");
rawohnepath = strcat(basepath,"mat_files/DataFileOfUnpreparedMeasurement_Volunteer2.mat");
rawmitpath = strcat(basepath,"mat_files/DataFileOfIRpreparedMeasurement_Volunteer2.mat");
trajpath = strcat(basepath,"../traj.mat");
baseoutbasepath = strcat(path_prefix, "/recons/volunteer2_");
rotatedeg = 131;

% elseif number == 3
% ...

end
