
% collects all results (M0,Mss,T1star,T1) from fit models (a),(b),(c) 

addpath(genpath("./lib/includes"))
addpath(genpath("./map_codes"))
addpath(genpath("./reco_scripts_for_paper"))

%%%%%%%%%%%%%%%%%%%% CONFIG
path_prefix = "";%/";
confnrs = 63;
%output_folder = "with results";
output_folder = "C://Temp/recons/";
slices = 1:35;
t1cutoffvalue = 3000; % cut off t1 values above this
%%%%%%%%%%%%%%%%%%%%


for confnr = confnrs

    [mask_path,rawohnepath,rawmitpath,trajpath,inbasepath,rotatedeg] = get_config(confnr, path_prefix);

    pathparts =  split(inbasepath, "/");
    namebeginning =  pathparts(end);

    if output_folder == "with results"
        outpath = strcat(inbasepath, "all_output.mat");
    else
        outpath = strcat(output_folder, namebeginning, "all_output.mat");

    end


    allout = zeros(size(slices,2),3,5,512,512);
    % dimensions: slices x fittype (combined, mit, ohne) x imtype (T1, T1star, Mss, M0_mit, M0_ohne) x ximsize x yimsize
    for sl = slices
        inpath = strcat(inbasepath, num2str(sl), "_end.mat");
        if ~isfile(inpath)
            if isfile(mask_path) % otherwise we are probably trying to generate masks right now
                load(mask_path)
                if sum(sum(squeeze(mask_all(sl,:,:)))) ~= 0
                    disp(strcat("MAP reconstructions of confnr ", num2str(confnr), " slice ",num2str(sl)," are not completed, using zeros as dummy"))
                end % else the reco was just skipped because mask said so
            end 
            continue
        end
        load(inpath)

        if allErgSep_comb == 0
            if isfile(mask_path) % otherwise we are probably trying to generate masks right now
                load(mask_path)
                if sum(sum(squeeze(mask_all(sl,:,:)))) ~= 0
                    disp(strcat("MAP reconstructions of confnr ", num2str(confnr), " slice ",num2str(sl)," are not completed, using zeros as dummy"))
                end % else the reco was just skipped because mask said so
            end 
            continue
        end

        %CombinedExpFit:
        Mss_c = squeeze(allErgSep_comb(1,1,:,:)); % = Mss_ohne
        M0_mit_c = squeeze(allErgSep_comb(1,2,:,:));
        M0_ohne_c = squeeze(allErgSep_comb(1,4,:,:) + allErgSep_comb(1,1,:,:));
        T1star_c = squeeze(allErgSep_comb(1,3,:,:)); % = T1star_ohne
        T1 = (T1star_c .* M0_ohne_c) ./ Mss_c ;
        T1(T1 > t1cutoffvalue) = t1cutoffvalue;

        allout(sl, 1, 1,:,:) = T1;
        allout(sl, 1, 2,:,:) = T1star_c;
        allout(sl, 1, 3,:,:) = Mss_c;
        allout(sl, 1, 4,:,:) = M0_mit_c;
        allout(sl, 1, 5,:,:) = M0_ohne_c;

        ll = size(allErgSep_uncomb,3);
        ur = size(allErgSep_uncomb,4);

        %UNCombinedExpFit:
        Mss_mit = zeros(512,512);
        Mss_mit(1:ll, 1:ur) =  squeeze(allErgSep_uncomb (:,1,:,:,2));
        M0_mit = zeros(512,512);
        M0_mit(1:ll, 1:ur) =  squeeze(allErgSep_uncomb(:,2,:,:,2));
        T1star_mit = zeros(512,512);
        T1star_mit(1:ll, 1:ur) = squeeze(allErgSep_uncomb(:,3,:,:,2));
        T1_mit = (T1star_mit .* M0_mit) ./ Mss_mit ;
        T1_mit(T1_mit > t1cutoffvalue) = t1cutoffvalue;


        allout(sl, 2, 1,:,:) = T1_mit;
        allout(sl, 2, 2,:,:) = T1star_mit;
        allout(sl, 2, 3,:,:) = Mss_mit;
        allout(sl, 2, 4,:,:) = M0_mit;
        allout(sl, 2, 5,:,:) = zeros(512,512);

        %UNCombinedExpFit:
        Mss_ohne = zeros(512,512);
        Mss_ohne(1:ll, 1:ur) =   squeeze(allErgSep_uncomb(:,1,:,:,1));
        M0_ohne = zeros(512,512);
        M0_ohne(1:ll, 1:ur) =  squeeze(allErgSep_uncomb(:,2,:,:,1) +  allErgSep_uncomb(:,1,:,:,1));
        T1star_ohne = zeros(512,512);
        T1star_ohne(1:ll, 1:ur) = squeeze(allErgSep_uncomb(:,3,:,:,1));
        T1_ohne = (T1star_ohne .* M0_ohne) ./ Mss_ohne ;
        T1_ohne(T1_ohne > t1cutoffvalue) = t1cutoffvalue;

        allout(sl, 3, 1,:,:) = T1_ohne;
        allout(sl, 3, 2,:,:) = T1star_ohne;
        allout(sl, 3, 3,:,:) = Mss_ohne;
        allout(sl, 3, 4,:,:) = zeros(512,512);
        allout(sl, 3, 5,:,:) = M0_ohne;

    end

    allrecos = zeros(size(slices,2),512,512);
    for sl = slices
        recopath = strcat(inbasepath, num2str(sl), "_comb0_reco.mat");
        if ~isfile(recopath)
            % disp(strcat("temporal average reconstructions of confnr ", num2str(confnr), " slice ",num2str(sl)," are not completed, using zeros as dummy"))
            continue
        end
        load(recopath);
        allrecos(sl,:,:) = temp_avg_reco;
    end

    %xtv(squeeze(allrecos(:,:,:)))

    save(outpath, 'allout', "allrecos");
end
% dimensions: slices x fittype (combined, mit, ohne) x imtype (T1, T1star, Mss, M0_mit, M0_ohne) x ximsize x yimsize

