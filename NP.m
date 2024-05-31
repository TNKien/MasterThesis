listPat         = [54, 57, 59,  61,   63,  69, 74, 77] ;
startSeiz       = [20 ,49 ,109 ,128 , 3  , 23, 47, 10] ;
durationSeiz    = [4  ,1  ,1   ,2   ,2   , 1 , 9 , 3  ] ;
startSeiz(startSeiz > 30) = 30 ;
startSeiz = startSeiz + 1 +1  ;
maxEdge = 0;
minEdge = 10000 ;
maxStr = 0;
minStr = 10000 ;
maxClust = 0;
minClust = 10000 ;

p = 77 ; 
indicator = 3;
bw = 6;
option = "" ;

index = find(listPat == p) ;
seiz = startSeiz(index) + 2  ;
pre = seiz - 5 ;
% post = seiz + durationSeiz(index) ;

% Création de la matrice 19x19 àpd d'une 18x18
%Pre
connMatrix_pre = ph(pre).Time ;               % load the 5D matrices (8x3x6x18x18)
sizeMatrix = size(connMatrix_pre) ; sizeMatrix(4:5) = 19 ;
newMat_pre = zeros(sizeMatrix);               % create the pre-seizure matrix of size 8x3x6x19x19
newMat_pre(:,:,:,1:7, 1:7)  = connMatrix_pre(:,:,:,1:7, 1:7); newMat_pre(:,:,:,1:7, 9:19) = connMatrix_pre(:,:,:,1:7, 8:18);
newMat_pre(:,:,:,9:19, 1:7) = connMatrix_pre(:,:,:,8:18, 1:7);newMat_pre(:,:,:,9:19,9:19) = connMatrix_pre(:,:,:,8:18,8:18);

coeff_pre = zeros(19,1) ;
coeff_pre(1:7) = ph_para(pre).Time(bw).bandWidth(indicator).Indicator{index,1}(1:7) ;
coeff_pre(9:19) = ph_para(pre).Time(bw).bandWidth(indicator).Indicator{index,1}(8:18) ;
str_pre = zeros(19,1) ;
str_pre(1:7) = ph_para(pre).Time(bw).bandWidth(indicator).Indicator{index,3}(1:7) ;
str_pre(9:19) = ph_para(pre).Time(bw).bandWidth(indicator).Indicator{index,3}(8:18) ;

aij_pre = squeeze(newMat_pre(index,indicator,bw,:,:)) ; ijw_pre = adj2edgeL(triu(aij_pre));
if max(aij_pre(aij_pre>0),[],'all') > maxEdge
    maxEdge = max(aij_pre(aij_pre>0),[],'all') ; end
if min(aij_pre(aij_pre>0),[],'all') < minEdge
    minEdge = min(aij_pre(aij_pre>0),[],'all') ; end
if min(coeff_pre(coeff_pre > 0)) < minClust
    minClust = min(coeff_pre(coeff_pre > 0)) ; end
if max(coeff_pre(coeff_pre > 0)) > maxClust
    maxClust = max(coeff_pre(coeff_pre > 0)) ; end
if min(str_pre(str_pre > 0)) < minStr
    minStr = min(str_pre(str_pre > 0)) ; end
if max(str_pre(str_pre > 0)) > maxStr
    maxStr = max(str_pre(str_pre > 0)) ; end

% %Post
% connMatrix_post = ph(post).Time ; newMat_post = zeros(sizeMatrix);               % create the post-seizure matrix of size 8x3x6x19x19
% newMat_post(:,:,:,1:7, 1:7)  = connMatrix_post(:,:,:,1:7, 1:7); newMat_post(:,:,:,1:7, 9:19) = connMatrix_post(:,:,:,1:7, 8:18);
% newMat_post(:,:,:,9:19, 1:7) = connMatrix_post(:,:,:,8:18, 1:7);newMat_post(:,:,:,9:19,9:19) = connMatrix_post(:,:,:,8:18,8:18);
% 
% coeff_post = zeros(19,1) ;
% coeff_post(1:7) = ph_para(post).Time(bw).bandWidth(indicator).Indicator{index,1}(1:7) ;
% coeff_post(9:19) = ph_para(post).Time(bw).bandWidth(indicator).Indicator{index,1}(8:18) ;
% str_post = zeros(19,1) ;
% str_post(1:7) = ph_para(post).Time(bw).bandWidth(indicator).Indicator{index,3}(1:7) ;
% str_post(9:19) = ph_para(post).Time(bw).bandWidth(indicator).Indicator{index,3}(8:18) ;
% 
% aij_post = squeeze(newMat_post(index,indicator,bw,:,:)) ; ijw_post = adj2edgeL(triu(aij_post));
% if max(aij_post(aij_post>0),[],'all') > maxEdge
%     maxEdge = max(aij_post(aij_post>0),[],'all') ; end
% if min(aij_post(aij_post>0),[],'all') < minEdge
%     minEdge = min(aij_post(aij_post>0),[],'all') ; end
% if min(coeff_post(coeff_post > 0)) < minClust
%     minClust = min(coeff_post(coeff_post > 0)) ; end
% if max(coeff_post(coeff_post > 0)) > maxClust
%     maxClust = max(coeff_post(coeff_post > 0)) ; end
% if min(str_post(str_post > 0)) < minStr
%     minStr = min(str_post(str_post > 0)) ; end
% if max(str_post(str_post > 0)) > maxStr
%     maxStr = max(str_post(str_post > 0)) ; end

%Seiz
connMatrix_seiz = ph(startSeiz(index)).Time ; newMat_seiz = zeros(sizeMatrix);% create the seizure matrix of size 8x3x6x19x19
newMat_seiz(:,:,:,1:7, 1:7)  = connMatrix_seiz(:,:,:,1:7, 1:7); newMat_seiz(:,:,:,1:7, 9:19) = connMatrix_seiz(:,:,:,1:7, 8:18);
newMat_seiz(:,:,:,9:19, 1:7) = connMatrix_seiz(:,:,:,8:18, 1:7);newMat_seiz(:,:,:,9:19,9:19) = connMatrix_seiz(:,:,:,8:18,8:18);

coeff_seiz = zeros(19,1) ;
coeff_seiz(1:7) = ph_para(seiz).Time(bw).bandWidth(indicator).Indicator{index,1}(1:7) ;
coeff_seiz(9:19) = ph_para(seiz).Time(bw).bandWidth(indicator).Indicator{index,1}(8:18) ;
str_seiz = zeros(19,1) ;
str_seiz(1:7) = ph_para(seiz).Time(bw).bandWidth(indicator).Indicator{index,3}(1:7) ;
str_seiz(9:19) = ph_para(seiz).Time(bw).bandWidth(indicator).Indicator{index,3}(8:18) ;

aij_seiz = squeeze(newMat_seiz(index,indicator,bw,:,:)) ; ijw_seiz = adj2edgeL(triu(aij_seiz));
if max(aij_seiz(aij_seiz>0),[],'all') > maxEdge
    maxEdge = max(aij_seiz(aij_seiz>0),[],'all') ; end
if min(aij_seiz(aij_seiz>0),[],'all') < minEdge
    minEdge = min(aij_seiz(aij_seiz>0),[],'all') ; end
if min(coeff_seiz(coeff_seiz > 0)) < minClust
    minClust = min(coeff_seiz(coeff_seiz > 0)) ; end
if max(coeff_seiz(coeff_seiz > 0)) > maxClust
    maxClust = max(coeff_seiz(coeff_seiz > 0)) ; end
if min(str_seiz(str_seiz > 0)) < minStr
    minStr = min(str_seiz(str_seiz > 0)) ; end
if max(str_seiz(str_seiz > 0)) > maxStr
    maxStr = max(str_seiz(str_seiz > 0)) ; end


% Plot the 3 networks
w_bounds = [minEdge, maxEdge] ;
% n_bounds = [minStr, maxStr] ;      %if we look at node Strenght
n_bounds = [minClust, maxClust] ;      %if we look at node Clustering coefficient
disp(w_bounds) ; 
disp(n_bounds) ;

%Pre
figure() ;
% subplot(1,2,1) ;
if option == "node" 
    f_PlotEEG_BrainNetwork(w_bounds,n_bounds, 19, ijw_pre, 'w_intact', coeff_pre, "n_nn2nx", "nocb"); 
else 
    f_PlotEEG_BrainNetwork(w_bounds,n_bounds, 19, ijw_pre, 'w_intact') ; 
end
titleNetwork = "5 minutes before seizure onset : patient n°" + num2str(p) ;
title(titleNetwork, FontSize=16) ;
%During
figure() ;
% subplot(1,2,2) ;
if option == "node"
    f_PlotEEG_BrainNetwork(w_bounds,n_bounds, 19, ijw_seiz, 'w_intact', coeff_seiz, "n_nn2nx", "nocb");
else 
    f_PlotEEG_BrainNetwork(w_bounds,n_bounds, 19, ijw_seiz, 'w_intact');
end
titleNetwork = "At seizure onset : patient n°" + num2str(p) ;
title(titleNetwork,FontSize=16) ;
%Post
% figure() ;
% if option == "node"
%     f_PlotEEG_BrainNetwork(w_bounds,n_bounds, 19, ijw_post, 'w_intact', coeff_post, "n_nn2nx", "nocb");
% else
%     f_PlotEEG_BrainNetwork(w_bounds,n_bounds, 19, ijw_post, 'w_intact') ;
% end
% titleNetwork = "Post-ictal : patient n°" + num2str(p) ;
% title(titleNetwork) ;