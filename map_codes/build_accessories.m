function [FT_op, mask_all, angle_all, cmap_all, wfull_all, temp_avg_reco] = build_accessories(raw, traj, w, mask_all,rotatedeg, onlyreco)
%BUILD_ACCESSORIES Builds all the necessary variables for the
%reconstruction from the raw data.

% expects raw as size 1 x coils x numspokes x spokelength
% mask has shape  1 x gridsize1 x gridsize2

incr = 1;
numspokes = size(raw,3);

%% Load/build trajectory

traj = traj(1:numspokes,:);
w = w(1:numspokes,:);

% include a rotation if necessary to ensure AP is up-down
if rotatedeg ~= 0
traj = traj.*exp(-1j*(rotatedeg/360)*2*pi ); %  rotate traj such that the head looks up
end

gridsize1=size(traj_frame,2)/4;
gridsize2=gridsize1;

%% Build FT operator %*GRIDDING*
disp("\n")
for frames = 1:incr:numspokes
    disp("Building FFT operator "+frames+"/"+numspokes)
    traj_frame = traj(frames:(frames+incr)-1,:);
    N = [size(traj_frame,2)/4,size(traj_frame,2)/4];
    ph = ones(size(traj_frame,2)/4);
    spokes = (frames:(frames+incr)-1);
    tr = traj(spokes,:);
    neww = w(frames:(frames+incr)-1,:);
    FT_op{(frames-1)/incr + 1} = BUILD THE GRIDDING OPERATOR *GRIDDING*
end
clear frames traj_frame N ph spokes tr neww

%% compute fully sampled temporal average for sensitivity maps, angle and
% wfull
coils = size(raw,2);
startspokeforfullreco = 250; % first line to be included in full reconstruction, should be after zero-crossing of all tissues

spokewise_reco = zeros(numspokes,coils,gridsize1,gridsize2); % it does not matter if this array is too big, we only add 0s to tempAv afterwards
% *GRIDDING*
for frames = startspokeforfullreco:incr:numspokes
    for coil = 1:coils
        spokewise_reco((frames-1)/incr + 1,coil,:,:) = FT_op{(frames-1)/incr + 1}'*squeeze(raw(1, coil, frames:(frames+incr)-1,:));
    end
    disp("Reconstructing raw data via nuFFT "+frames+"/"+numspokes)
end
tempAv = squeeze(sum(spokewise_reco,1));

%% ADD: coil sensitivities/ function to determine coil sensitivities. *COIL COMBINATION*
% e.g. [rec,cmap,wfull] = *COIL COMBINATION*(tempAv);
% expected output/data
% rec = coil combined image with phase information, size(rec)= [gridsize1 gridsize2];
% cmap = coil maps size(cmap)=[nb of coils gridsize1, gridsize2];


angle_all = angle(rec);
cmap_all = cmap;
wfull_all = wfull;

temp_avg_reco = rec;

if onlyreco
    return
end

%%

% insert leading singleton dimensions
angle_all = reshape(angle_all, [1,size(angle_all)]);
cmap_all = reshape(cmap_all, [1,size(cmap_all)]);
wfull_all = reshape(wfull_all, [1,size(wfull_all)]);

end