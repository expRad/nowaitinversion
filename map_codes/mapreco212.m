function [allErgSep] = mapreco212(rawOhne, rawMit, combined, outbasepath, traj, w, sl_mask_all, blocknr,rotatedeg, onlyreco)
% expects rawOhne and rawMit as 1 x coils x numspokes x spokelength

[FT_op, mask_all, angle_all, cmap_all, wfull_all, temp_avg_reco] = build_accessories(rawOhne, traj, w, sl_mask_all,rotatedeg, onlyreco);

if outbasepath ~= "" % because in this case, we assume that no intermediate steps should be saved
    outbasepath = strcat(outbasepath, "comb", num2str(combined), "_");
    save(strcat(outbasepath,"reco.mat"),"temp_avg_reco")
end

if onlyreco
    allErgSep = 0;
    return
end

tic
% use dc after zero as initialization
ite=0;

% this is just fixed
%% ADD you sequence settings here
% TR - repetition time
% t0 - time delay between start of inversion recovery (e.g. centre of
% inversion pulse) and start of unprepared transition phase (e.g. first
% slice excitation pulse)
incr = 1;
TR = 7.5*incr;
t0 = TR + 23;
numspokes = size(rawOhne,3);

z = zeros(1,numspokes,512,512);
[dc_mit_im2,~] = dataConsis_MitIR_SingSlc(rawMit,FT_op,incr,numspokes,z,cmap_all,wfull_all,ite);
[dc_ohne_im2,~] = dataConsis_OhneIR_SingSlc(rawOhne,FT_op,incr,numspokes,z,cmap_all,wfull_all,ite);

savepath = ""; % uncomment this to save only the result of the last step
%savepath = outbasepath;  % uncomment this to save all intermediate steps


for i = 1:3:blocknr*3
    ite = i
    
    [fit_mit_im1,fit_ohne_im1, dc_mit_im1, dc_ohne_im1, ~] = map_step(dc_mit_im2,dc_ohne_im2, ite, savepath, mask_all, TR, t0, angle_all, rawMit, rawOhne, FT_op, incr, numspokes, cmap_all, wfull_all, combined,1);
    
    ite = i+1
    
    [fit_mit_i,fit_ohne_i, dc_mit_i, dc_ohne_i, ~] = map_step(dc_mit_im1,dc_ohne_im1, ite, savepath, mask_all, TR, t0, angle_all, rawMit, rawOhne, FT_op, incr, numspokes, cmap_all, wfull_all, combined,1);
    
    ite = i+2
    
    [fit_mit_ip1,fit_ohne_ip1,dc_mit_ip1,dc_ohne_ip1] = linear_step_framewise(ite,fit_mit_im1, fit_mit_i,fit_ohne_im1,fit_ohne_i,dc_mit_im1,dc_mit_i,dc_ohne_im1,dc_ohne_i, savepath);
    

    dc_mit_im2 = dc_mit_ip1;
    dc_ohne_im2 = dc_ohne_ip1;
    clear fit_mit_i fit_ohne_i dc_mit_i dc_ohne_i fit_mit_im1 fit_ohne_im1 dc_mit_im1 dc_ohne_im1
    
end

ite = blocknr*3+1

%%[fit_mit_im1,fit_ohne_im1, dc_mit_im1, dc_ohne_im1, ~] = map_step(dc_mit_im2,dc_ohne_im2, ite, savepath, mask_all, TR, t0, angle_all, rawMit, rawOhne, FT_op, incr, numspokes, cmap_all, wfull_all, combined,1);
[fit_mit_im1,fit_ohne_im1, dc_mit_im1, dc_ohne_im1, ~] = map_step(dc_mit_im2,dc_ohne_im2, ite, outbasepath, mask_all, TR, t0, angle_all, rawMit, rawOhne, FT_op, incr, numspokes, cmap_all, wfull_all, combined,1);

ite = blocknr*3+2

[fit_mit_i,fit_ohne_i, dc_mit_i, dc_ohne_i, allErgSep] = map_step(dc_mit_im1,dc_ohne_im1, ite, outbasepath, mask_all, TR, t0, angle_all, rawMit, rawOhne, FT_op, incr, numspokes, cmap_all, wfull_all, combined,0);

toc
end