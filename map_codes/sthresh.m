
 
function filim=sthresh(orgim,thresh)
 
       
       X=orgim;
       T=thresh;
       
%     ind=find(abs(X)<=T);
%     ind1=find(abs(X)>T);
    X(find(abs(X)<=T))=0;
    X(find(abs(X)>T))=sign(X(find(abs(X)>T))).*(abs(X(find(abs(X)>T)))-T);
    filim=X;
end
 