clear; clc;

dataFolder = pwd + "\Data crises\" ;
nbPatients  = [ 54    57    59    61    63    69    74    77];
start       = [ 20,   49,   109,  128,  3,    23,   47,   314]*250*60;       % seizure start time for each patient
nbWindow    = 60 ;                                                           % # of sliding windows
duration    = 15000 ;                                                        % duration of sliding windows (1min)

channels    = ["F3","C3","P3","Cz","F4","C4","P4","Fp1","Fp2","F7","T3","T5","O1","O2","T6","T4","F8","Fz","Pz"];
ref         = "C3";
reref       = "";
fhband      = [1,4,8,12,1,1];      % delta, theta, alpha, beta, quasi-broad, broad
flband      = [4,8,12,20,12,20];
window      = 1;                   % in seconds
overlap     = 0.5;                 % in seconds 

corr = struct() ;
ph = struct() ;

for i = 1:nbWindow
    corr(i).Time = zeros(length(nbPatients),2,length(fhband),length(channels)-1,length(channels)-1);    % un tableau 5D pour correlation
    ph(i).Time = zeros(length(nbPatients),3,length(fhband),length(channels)-1,length(channels)-1);      % un tableau 5D pour la phase
end 

corr_para = struct() ;
ph_para = struct() ;

for p = nbPatients
    disp("Processing patient n° " + p)
    load(dataFolder + p + ".mat")

    index = find(nbPatients == p) ;

    t = max(1, start(index)+ 1 -30*250*60) ;  % the first sliding window start 30min before the seizure (if possible)

    
        for i = 1:nbWindow
            if sum(sum(EEG.data(1).series(t:t+duration-1 , :))) ~= 0
                disp("Interval n° " + i)
                disp(t) ;
                
                corr(i).Time(index,:,:,:,:) = connectivity(EEG,t,duration,channels,ref,reref,fhband,flband,window,overlap,correlation(50));
                ph(i).Time(index,:,:,:,:) = connectivity(EEG,t,duration,channels,ref,reref,fhband,flband,window,overlap,phase);
                
                corr_para(i).Minute = i ; ph_para(i).Minute = i ;
                
                for b = 1:length(fhband)
                    for ind = 1:3
                        if ind ~= 3
                            corr_para(i).Time(b).bandWidth(ind).Indicator(index,:) = parameters(corr(i).Time(index,ind,b,:,:)) ; end
                        ph_para(i).Time(b).bandWidth(ind).Indicator(index,:) = parameters(ph(i).Time(index,ind,b,:,:)) ;
                        
                    end
                    
                end
            end
            t = t + duration ;
        end
   
end

save(pwd + "\Results\corr.mat", "corr");
save(pwd + "\Results\ph.mat", "ph");
save(pwd + "\Results\corr_para.mat", "corr_para");
save(pwd + "\Results\ph_para.mat", "ph_para");