function [dataConsisOhne,residualOhne] = dataConsis_OhneIR_SingSlc(raw,FT_op,incr,nArms,fit,cmap,wfull,ite)   

    % This setup works for reconstrctions schemes with ONE spiral arm
    % per image (per frame)

    % Attention should be given to which dataset is the input (Mit or Ohne
    % IR) to adjust the name of the saving file accordingly)

    residualOhne = zeros(size(fit,1),size(fit,2),size(fit,3),size(fit,3));
    dataConsisOhne = zeros(size(fit,1),size(fit,2),size(fit,3),size(fit,3));
    
for slice = 1:size(raw,1)
  
    disp(slice)

    multiCoil_resid = zeros(1,size(cmap,2),size(fit,3),size(fit,3));
    resid = zeros(size(fit,2),size(fit,3),size(fit,3));
    
    newfit = squeeze(fit(slice,:,:,:));
    newcmap = squeeze(cmap(slice,:,:,:));
    newwfull = squeeze(wfull(slice,:,:,:));
    newraw = squeeze(raw(slice,:,:,:));

    for frames = 1:incr:nArms
        
        disp(strcat("DataConsis_OhneIR ite=", num2str(ite)," ", num2str(frames), " / ", num2str(nArms)))

        undersamp = newraw(:,frames:(frames+incr)-1,:);

        for coils = 1:size(raw,2)

                sc_frame = newcmap(coils,:,:).*newfit((frames-1)/incr+1,:,:);

                k_mod(coils,:,:) = FT_op{(frames-1)/incr + 1} *squeeze(sc_frame); % *GRIDDING*
                    
                multiCoil_resid(:,coils,:,:) = FT_op{(frames-1)/incr + 1}'*squeeze(k_mod(coils,:,:)  - undersamp(coils,:,:)); % *GRIDDING*
               
        end           
          
               resid((frames-1)/incr + 1,:,:) = squeeze((sum(newwfull.* squeeze(multiCoil_resid(1,:,:,:)),1)));

   
    end
 
    result_dtCons = newfit - resid;
    
    residualOhne(slice,:,:,:) = resid;
    dataConsisOhne(slice,:,:,:) = result_dtCons;
    
    clear resid result_dtCons newfit newraw newcmap newwfull k_mod sc_frame undersamp multiCoil_resid
end

end
  