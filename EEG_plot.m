a = 1 ; %parfois, certaines data ont 2 tableaux jsp pq
fs = 250 ; 
listPat         = [54 ,57 ,59  ,61  , 63 , 69, 74, 77 ] ;
startSeiz       = [20 ,49 ,109 ,128 , 3  , 23, 47, 10] ;
durationSeiz    = [4  ,1  ,1   ,2   , 2   , 1 , 9 , 3  ] ;

p = 8 ;

sizeData = size(EEG.data(a).series) ;
temps = 0:60*60*250 - 1 ;
temps = temps/250/60;
start = (max(0,startSeiz(p)-30)*60*250 ) + 1 ; endTime = start + 60*60*250 - 1 ;
figure ;
for j = 1
    % subplot(7,3,j) ;
    plot(temps, EEG.data(a).series(start:endTime,j)) ;
    hold on ;
    % rectangle('Position', [min(startSeiz(p),30), min(EEG.data(a).series(start:endTime,j)), durationSeiz(p),max(EEG.data(a).series(start:endTime,j))-min(EEG.data(a).series(start:endTime,j)) ], 'FaceColor',[0.5 1 0.5 0.5], 'EdgeColor', 'none') ;
    rectangle('Position', [min(startSeiz(p),30), -500, durationSeiz(p),1000 ], 'FaceColor',[0.5 1 0.5 0.5], 'EdgeColor', 'none') ;
    xline(10, 'r--', 'LineWidth', 1.5); 
    xline(5, 'r--', 'LineWidth', 1.5); 
    hold off ;
    Title = "Electrode " + EEG.channels(j).Name ;
    xlabel("Time (min)", FontSize=16) ;
    ylabel("Voltage (µV)", FontSize=16) ;
    title(Title, FontSize=16) ; 
    % xlabel("Time (min)") ;
    % ylabel("Voltage (µV)") ;
    % title(Title) ;
    ylim([-500, 500]) ;
    % axis tight ;
end