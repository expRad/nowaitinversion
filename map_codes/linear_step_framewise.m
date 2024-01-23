function [fit_mit_ip1,fit_ohne_ip1,dc_mit_ip1,dc_ohne_ip1] = linear_step_framewise(ite,fit_mit_im1, fit_mit_i,fit_ohne_im1,fit_ohne_i,dc_mit_im1,dc_mit_i,dc_ohne_im1,dc_ohne_i,outbasepath)
fit_mit_i = squeeze(fit_mit_i);
fit_ohne_i = squeeze(fit_ohne_i);
dc_mit_i = squeeze(dc_mit_i);
dc_ohne_i = squeeze(dc_ohne_i);

fit_mit_im1 = squeeze(fit_mit_im1);
fit_ohne_im1 = squeeze(fit_ohne_im1);
dc_mit_im1 = squeeze(dc_mit_im1);
dc_ohne_im1 = squeeze(dc_ohne_im1);

% linear step
[fit_mit_ip1, dc_mit_ip1, ~, ~] = line_step_framewise(fit_mit_im1, fit_mit_i-fit_mit_im1, dc_mit_im1, dc_mit_i-dc_mit_im1, 1, true);
[fit_ohne_ip1, dc_ohne_ip1, ~, ~] = line_step_framewise(fit_ohne_im1, fit_ohne_i-fit_ohne_im1, dc_ohne_im1, dc_ohne_i-dc_ohne_im1, 1, true);

if outbasepath ~= ""
fit_ohne = squeeze(fit_ohne_ip1);
fit_mit = squeeze(fit_mit_ip1);
dc_mit = squeeze(dc_mit_ip1);
dc_ohne = squeeze(dc_ohne_ip1);
save(strcat(outbasepath, num2str(ite), ".mat"), "fit_ohne", "fit_mit", "dc_mit", "dc_ohne");
end
end

