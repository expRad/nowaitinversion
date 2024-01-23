addpath(genpath("./lib/includes"))
addpath(genpath("./map_codes"))
addpath(genpath("./reco_scripts_for_paper"))

%%%%%%%%%%%%%%%%%%%% CONFIG
numspokes = 400;  % spiral spokes of the measurement. 
blocknr = 2;      % how many blocks of 2 map steps followed by 1 linear step should be done. we decided on 2 for the paper.
onlyreco = false;     % set to 1 if you want to do only temporal average recos, and 0 if you want to do MAP reconstruction for quantitative maps
onlycomb0 = false;    % set to 1 if you want to run only the combined reconstruction (model (c)) and 0 for all reconstructions (models (a),(b) and (c))
confnrs = [1 2]; % which datasets to run. dataset-specific configs are set in get_config.m
slices = 1:35;     % which slices to reconstruct. can be any subset of all slices in any order.
path_prefix = ""; % this is prefixed to all paths for loading and saving data. 
%%%%%%%%%%%%%%%%%%%

for confnr = confnrs 

    [mask_path,rawohnepath,rawmitpath,trajpath,baseoutbasepath,rotatedeg] = get_config(confnr, path_prefix);

    load(rawohnepath)
    rawOhne_allslc = raw(:,:,1:numspokes,:);
    clear raw
    load(rawmitpath);
    rawMit_allslc = raw(:,:,1:numspokes,:);
    clear raw
    allslices = size(rawMit_allslc,1);
    if ~onlyreco
        load(mask_path);
        if confnr == 7
            mask_all = maskOhne;
        end
    else
        mask_all = ones(allslices,1,1);
    end
    load(trajpath);

    for slice=slices % 1:allslices

        outbasepath = strcat(baseoutbasepath, num2str(slice),"_");
        sl_mask_all = mask_all(slice,:,:);

        if ~(max(max(sl_mask_all)) == 0)
            rawOhne = rawOhne_allslc(slice,:,:,:);
            rawMit = rawMit_allslc(slice,:,:,:);

            %%% Run reconstruction
            [allErgSep_uncomb] = mapreco212(rawOhne, rawMit, 0, outbasepath, traj, w, sl_mask_all, blocknr,rotatedeg,onlyreco); % Uncombined fit (Methods (a) and (b) )
            if ~onlycomb0
                [allErgSep_comb] = mapreco212(rawOhne, rawMit, 1, outbasepath, traj, w, sl_mask_all, blocknr,rotatedeg,onlyreco); % Combined fit (Method (c) )
            else
               allErgSep_comb = 0;
            end
            save(strcat(outbasepath, "end", ".mat"), "allErgSep_uncomb", "allErgSep_comb");

        end

        disp(slice)

    end

end


