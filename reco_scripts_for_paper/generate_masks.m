%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIFIC MASKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_prefix = "";
confnr = 63;

[mask_path,rawohnepath,rawmitpath,trajpath,baseoutbasepath,rotatedeg] = get_config(confnr, path_prefix);
recopath = strcat(baseoutbasepath, "all_output.mat");
load(recopath)
mask_all = zeros(size(allrecos,1),512,512);
for sl=1:size(allrecos,1)
    image = squeeze(abs(allrecos(sl,:,:)));
    Im1=n2one(image);
    Im1(isnan(Im1)) = 1; 
    Mask = logical(sthresh( Im1 , 0.2)); % adjust this value as needed
    Mask = imfill(Mask,'holes');
    if max(max(abs(allrecos(sl,:,:)))) < 3e-4 % threshold, when to considers a slice empty; adjust this value as needed
        mask_all(sl,:,:) = 0;
    else
        mask_all(sl,:,:) = Mask;
    end
end

% might be worth to insert some kind of display, e.g. imagesc

save(mask_path, "mask_all")


