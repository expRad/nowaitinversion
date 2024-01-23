function [f3,d3, lamb, mu] = line_step_framewise(f1, f2, d1, d2, damping, compl)
% performs the linear step for each frame separately

% assumes size of f and d is
% numframes x ximsize x yimsize


[numframes, ximsize, yimsize] = size(f1);


f3 = zeros(numframes, ximsize, yimsize);
d3 = zeros(numframes, ximsize, yimsize);

for n=1:numframes
    [fn, dn,lamb,~] = line_step(squeeze(f1(n,:,:)),squeeze(f2(n,:,:)),squeeze(d1(n,:,:)),squeeze(d2(n,:,:)),damping, compl);
    f3(n,:,:) = fn;
    d3(n,:,:) = dn;
    disp(n)
    disp(lamb)
end
lamb = 1;
mu = 1;
end