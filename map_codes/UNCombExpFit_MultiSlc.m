
function [allErgSep_new_ohne,allErgSep_new_mit,fitMit,fitOhne] = UNCombExpFit_MultiSlc(s1Mit,s1Ohne,mask1,TR,t0,anglefit,ite)

% A(1) = Msteadystate
% A(2) = M0
% A(3) = T1stern
% A(4) = M0 - Mss
%

% fun1 - LL inversion recovery
fun1 = @(A,xdata) A(1)- (A(2)+A(1)).*exp(-xdata./A(3));

% fun3 - relaxation, no inversion, different M0 than with inversion % (imperfect inversion)

fun3 = @(A,xdata) A(1)+ (A(2)).*(exp(-xdata./A(3))); % A(2) stands for M0-Mss % changed for UNCombined

funSep = @(A,xdata) [fun1(A,xdata(:,1)),  fun3(A,xdata(:,2))]; % dfferent M0, A(1)-A(4)

nb = size(s1Mit,2);

xachs = t0:TR:((TR*nb)+t0)-1;
x_mtx = [xachs(:) xachs(:)-23]; % 23.6530 xachsMit xAchsOhne -shifted as relaxation starts after inversion or after first GRE-pulse

dim = size(s1Mit,3);

% init results from funSep

fitMit = zeros(size(s1Mit));
fitOhne = zeros(size(s1Ohne));
%allErgSep_new = zeros(24,4,dim,dim);

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

    allErgSep_mit = zeros(4,dim,dim);
    allErgSep_ohne = zeros(4,dim,dim);
    allErgSep_new_mit = zeros(size(s1Mit,1), 4, dim,dim);
    allErgSep_new_ohne = zeros(size(s1Mit,1), 4, dim,dim);

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

                UBSep = [5*M0_start 50*M0_start 7000 2000]; % upper bound [Mss M0 T1 M0_2]

                ErgSep_mit = lsqcurvefit(fun1, x0Sep, x_mtx(:,1), real(y_mtxa(:,1)),LBSep,UBSep,opts);
                ErgSep_ohne = lsqcurvefit(fun3, x0Sep, x_mtx(:,2), real(y_mtxa(:,2)),LBSep,UBSep,opts); % changed for UNCombined

                           ErgSep_curves_mit = fun1(ErgSep_mit,x_mtx(:,1));
                           ErgSep_curves_ohne = fun3(ErgSep_ohne,x_mtx(:,2)); % changed for UNCombined

                allErgSep_mit(:,la,lb) = ErgSep_mit;
                allErgSep_ohne(:,la,lb) = ErgSep_ohne;

                           s1MitSep_fit(:,la,lb) = ErgSep_curves_mit.*1e-6; % changed for UNCombined
                           s1OhneSep_fit(:,la,lb) = ErgSep_curves_ohne.*1e-6;  % undo scaling from before

                new_s1MitSep_fit(:,la,lb) = s1MitSep_fit(:,la,lb).*exp(1i*squeeze(nangle1(1,la,lb)));
                new_s1OhneSep_fit(:,la,lb) = s1OhneSep_fit(:,la,lb).*exp(1i*squeeze(nangle1(1,la,lb)));

            end

        end


    end

    clear MssMit_fit M0Ohne_fit M0Mit_fit T1startSep_fit clear s1MitSep_fit s1OhneSep_fit

    fitMit(slice,:,:,:) = new_s1MitSep_fit; clear new_s1MitSep_fit
    fitOhne(slice,:,:,:) = new_s1OhneSep_fit; clear new_s1OhneSep_fit

    allErgSep_new_mit(slice,:,:,:) = allErgSep_mit; clear allErgSep
    allErgSep_new_ohne(slice,:,:,:) = allErgSep_ohne; clear allErgSep


end

end

