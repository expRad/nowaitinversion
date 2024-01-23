function [fit_mit_i,fit_ohne_i, dc_mit_i, dc_ohne_i, allErgSep] = map_step(dc_mit_im1,dc_ohne_im1, ite, outbasepath, mask_all, TR, t0, angle_all, rawMit, rawOhne, FT_op, incr, numspokes, cmap_all, wfull_all, combined, do_dc)
if combined
[allErgSep,fit_mit_i,fit_ohne_i] = CombExpFit_MultiSlc(reshape(dc_mit_im1, [1,numspokes,512,512]),reshape(dc_ohne_im1, [1,numspokes,512,512]),mask_all,TR,t0,angle_all,ite);
else
[allErgSep_new_ohne,allErgSep_new_mit,fit_mit_i,fit_ohne_i] = UNCombExpFit_MultiSlc(reshape(dc_mit_im1, [1,numspokes,512,512]),reshape(dc_ohne_im1, [1,numspokes,512,512]),mask_all,TR,t0,angle_all,ite);
allErgSep = cat(5, allErgSep_new_ohne, allErgSep_new_mit);
end

if do_dc
[dc_mit_i,~] = dataConsis_MitIR_SingSlc(rawMit,FT_op,incr,numspokes,fit_mit_i,cmap_all,wfull_all,ite);
[dc_ohne_i,~] = dataConsis_OhneIR_SingSlc(rawOhne,FT_op,incr,numspokes,fit_ohne_i,cmap_all,wfull_all,ite);
else
    dc_mit_i=0;
    dc_ohne_i = 0;
end
   
if outbasepath ~= ""
fit_ohne = squeeze(fit_ohne_i);
fit_mit = squeeze(fit_mit_i);
dc_mit = squeeze(dc_mit_i);
dc_ohne = squeeze(dc_ohne_i);
save(strcat(outbasepath, num2str(ite), ".mat"), "fit_ohne", "fit_mit", "dc_mit", "dc_ohne", "allErgSep");
end
end

