% Global data
listPara = ["Global clustering coefficient", "Characteristic path length","Average node strength", "Global efficiency"] ;
listPat         = [54 ,57 ,59  ,61  , 63 , 69, 74, 77 ] ;
startSeiz       = [20 ,49 ,109 ,128 , 3  , 23, 47, 10] ;
durationSeiz    = [4  ,1  ,1   ,2   ,2   , 1 , 9 , 3  ] ;
numPara = length(listPara);
numPat = length(listPat) ;
bands = ["delta", "theta", "alpha", "beta", "quasi-broad", "broad"] ;
windowSize = 1 ;            % size of sliding window (in minutes)
time = 0:windowSize:59 ;
pval = zeros(6,4) ; pval_CCC = zeros(6,4) ;


for bw = 1
    % Load the parameters
    para_phase = zeros(size(ph_para,1),numPat,numPara) ; % create a 3D matrix 60x8x6 (nbInterval x nbPatient x nbPara)
    para_corr  = zeros(size(ph_para,1),numPat,numPara) ;
    for t = 1:size(ph_para,2)
        for p = 1:numPat
            for para = 1:numPara
                para_phase(t,p,para) = mean(ph_para(t).Time(bw).bandWidth(3).Indicator{p,para}) ;        % we focus only on wPLI for the phase measures
                para_corr(t,p,para) = mean(corr_para(t).Time(bw).bandWidth(2).Indicator{p,para});        % we focus only on CCC for the correlation measures
            end
        end
    end
    % Correction for patient 57 at 10th minute
    para_phase(10,2,:) = (para_phase(9,2,:) + para_phase(11,2,:)) / 2 ;
    para_corr(10,2,:) = (para_corr(9,2,:) + para_corr(11,2,:)) / 2 ;

    for p = 8
        figure() ;
        titleFig = "Patient n°" + listPat(p) + " : Global metrics evolution"  ;
        sgtitle(titleFig, Fontsize = 16) ;
        for j = 1:numPara
            subplot(2,2,j) ;
            % plot(time, para_phase(:,p,j), 'b') ;
            % hold on ;
            plot(time, para_corr(:,p,j), 'r');
            hold on ;

            % Computation of p-value
            start = min(30, startSeiz(p)) ;
            start = start + 1 ;
            data_outside_crisis = para_phase([1:start + 1,start + durationSeiz(p) + 1 :end],p,j) ;
            data_during_crisis = para_phase(start+1:start + durationSeiz(p),p,j);
            data_outside_crisis2 = para_corr([1:start + 1,start + durationSeiz(p) + 1 :end],p,j) ;
            data_during_crisis2 = para_corr(start+1:start + durationSeiz(p),p,j);
            % Call the permutation test
            [observed_diff, p_value] = stat_test(data_outside_crisis, data_during_crisis);
            [observed_diff2, p_value2] = stat_test(data_outside_crisis2, data_during_crisis2);
            pval(bw,j) = p_value ; pval_CCC(bw,j) = p_value2 ;


            min_ph = min(para_phase(:,p,j)) ; max_ph = max(para_phase(:,p,j));
            min_corr = min(para_corr(:,p,j)) ; max_corr = max(para_corr(:,p,j)) ;
            % rectangle('Position', [min(startSeiz(p),30), min(min_ph,min_corr), durationSeiz(p),max(max_ph,max_corr)-min(min_ph,min_corr) ], 'FaceColor',[0.5 1 0.5 0.5], 'EdgeColor', 'none') ;
            % rectangle('Position', [min(startSeiz(p),30), min(para_phase(:,p,j)), durationSeiz(p),max(para_phase(:,p,j))-min(para_phase(:,p,j)) ], 'FaceColor',[0.5 1 0.5 0.5], 'EdgeColor', 'none') ;
            rectangle('Position', [min(startSeiz(p),30), min(para_corr(:,p,j)), durationSeiz(p),max(para_corr(:,p,j))-min(para_corr(:,p,j)) ], 'FaceColor',[0.5 1 0.5 0.5], 'EdgeColor', 'none') ;
            hold off ;


            xlabel('Time (min)');
            ylabel('Parameter Value');
            titlePlot = listPara(j) ;
            %titlePlot = " p-value wPLI : " + p_value + " and p-value CCC : " +p_value2 ;
            title(titlePlot, FontSize=14);
            % ylim([min(para_corr(11:end,p,j)), max(para_corr(11:end,p,j))]) ;
            axis tight ;
            
        end
        pathName = "D:\Mémoire - codes finaux\Time-series\" + listPat(p) + "_timeseries_" + bands(bw) + ".jpg";
        figureHandle = gcf;
        set(figureHandle, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
        % saveas(figureHandle, pathName);
        break
    end
end
