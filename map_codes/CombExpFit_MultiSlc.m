function [allErgSep_new,fitMit,fitOhne] = CombExpFit_MultiSlc(s1Mit,s1Ohne,mask1,TR,t0,anglefit,ite)


% A(1) = Msteadystate
% A(2) = M0
% A(3) = T1star
% A(4) = M0 - Mss
%

% fun1 - LL inversion recovery
fun1 = @(A,xdata) A(1)- (A(2)+A(1)).*exp(-xdata./A(3));


% fun3 - relaxation, no inversion, different M0 than with inversion
% (imperfect inversion)
fun3 = @(A,xdata) A(1)+ (A(4)).*(exp(-xdata./A(3))); % A(4) stands for M0-Mss

funSep = @(A,xdata) [fun1(A,xdata(:,1)),  fun3(A,xdata(:,2))]; % dfferent M0, A(1)-A(4)

nb = size(s1Mit,2);

xachs = t0:TR:((TR*nb)+t0)-1;
x_mtx = [xachs(:) xachs(:)-23]; % 23.6530 xachsMit xAchsOhne -shifted as relaxation starts after inversion or after first GRE-pulse

dim = size(s1Mit,3);

% init results from funSep

fitMit = zeros(size(s1Mit));
fitOhne = zeros(size(s1Ohne));
allErgSep_new = zeros(size(s1Mit,1),4,dim,dim);

opts = optimset('Display','off');

%% do pixelwise fitting, slice by slice
for slice = 1:size(s1Mit,1)

       MssMit_fit = zeros(dim,dim); % A(1)
       M0Ohne_fit = zeros(dim,dim); % A(2)
       M0Mit_fit = zeros(dim,dim); %A(4)
       T1starSep_fit = zeros(dim,dim); % A(3)
    
    allErgSep = zeros(4,dim,dim);

    s1OhneSep_fit = zeros(nb,dim,dim); %fitted exp-curves
    s1MitSep_fit = zeros(nb,dim,dim); %fitted exp-curves
    new_s1MitSep_fit = zeros(nb,dim,dim); %fitted exp-curves with the angle
    new_s1OhneSep_fit = zeros(nb,dim,dim); %fitted exp-curves with the angle

    nipfit_Mit = squeeze(s1Mit(slice,:,:,:));
    nipfit_Ohne = squeeze(s1Ohne(slice,:,:,:));

    nangle1 = anglefit(slice,:,:);

    maskfit = mask1(slice,:,:);

    for la=1:dim

        disp(strcat("Fit ite=", num2str(ite)," ", num2str(la), " / ", num2str(dim)))

        for lb=1:dim

            if maskfit(1,la,lb)

                % signal for 1 px at position la,lb from both measurements
                y_mtx = [squeeze(nipfit_Mit(:,la,lb)) squeeze(nipfit_Ohne(:,la,lb))];

                y_mtx = y_mtx.*1e6; % this is a simple empirical upscaling, as matlab fitting sometimes struggles with small numbers 

                %rotate onto real axis
                y_mtxa = y_mtx.*exp(-1i*squeeze(nangle1(1,la,lb)));
 
                % initialize fitting parameters
                M0_start = abs(y_mtxa(1,2)); % first value of measuement without inversion
                MssStart = mean(mean(real(y_mtxa(floor(3*end/4):end,:)))); % mean of last values along relaxation curve
                MssStart = min([MssStart M0_start]); % change in case of bad noise in the end of short sigs or incomplete relaxation
                T1start = 1000;

                x0Sep = [MssStart M0_start T1start M0_start-MssStart]; % start values [Mss M0 T1 M0_2]

                LBSep = [0 0 50 0]; % lower bound [Mss M0 T1 M0_2]

                UBSep = [5*M0_start 50*M0_start 7000 50*M0_start];

                % perform fitting
                ErgSep = lsqcurvefit(funSep, x0Sep, x_mtx, real(y_mtxa),LBSep,UBSep,opts);

                         ErgSep_curves = funSep(ErgSep,x_mtx);

                         MssMit_fit(la,lb) = ErgSep(1);
                         M0Ohne_fit(la,lb) = ErgSep(2);
                         M0Mit_fit(la,lb) = ErgSep(4)+ErgSep(1);
                         T1starSep_fit(la,lb) = ErgSep(3);

                allErgSep(:,la,lb) = ErgSep;

                           s1MitSep_fit(:,la,lb) = ErgSep_curves(:,1).*1e-6; % undo scaling from before
                           s1OhneSep_fit(:,la,lb) = ErgSep_curves(:,2).*1e-6;

                new_s1MitSep_fit(:,la,lb) = s1MitSep_fit(:,la,lb).*exp(1i*squeeze(nangle1(1,la,lb)));
                new_s1OhneSep_fit(:,la,lb) = s1OhneSep_fit(:,la,lb).*exp(1i*squeeze(nangle1(1,la,lb)));

            end

        end


    end

    clear MssMit_fit M0Ohne_fit M0Mit_fit T1startSep_fit
    clear s1MitSep_fit s1OhneSep_fit

    %collect results from all px and slices
    fitMit(slice,:,:,:) = new_s1MitSep_fit; clear new_s1MitSep_fit
    fitOhne(slice,:,:,:) = new_s1OhneSep_fit; clear new_s1OhneSep_fit

    allErgSep_new(slice,:,:,:) = allErgSep; clear allErgSep
end

end

