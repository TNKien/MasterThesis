% Class Name:   data
% Description:  Data fetching and pre-processing
% Author:       Lise Cottin
% Date:         September 25, 2023

classdef data
    properties
        samples         {mustBeNumeric}
        channels
        fs              {mustBeNumeric}
        rawData         % EEG data [samples x channels]
        mastoids        % EEG mastoids data
        laplace = {["F3","Fp1","F7","Fz"],["P3","T5","O1","Pz"],["F4","Fp2","Fz","C4","F8"],["T5","T3","P3","O1"],...
                ["Cz","Fz","Pz","C4","F3","F4","P4","P3"],["C4","F4","Cz","P4","T4"],["P4","C4","Pz","O2","T6"],...
                ["Fp1","F7","F3","Fz","Fp2"],["F7","Fp1","F3","T3"],["Fp2","Fz","F4","F8","Fp1"],["T3","F7","T5","A1"],...
                ["O1","T5","P3","Pz","O2"],["O2","T6","P4","Pz","O1"],["T6","T4","P4","O2","C4"],["T4","F8","C4","T6","A2"],...
                ["F8","Fp2","F4","T4","C4"],["Fz","Fp1","Fp2","F3","Cz","F4"],["Pz","O1","O2","P3","Cz","P4"]};
    end

    methods
        function obj = data(EEG,start,samples,channels,ref)
            %if strcmp(samples,'all'); obj.samples = size(EEG.data,1);  %pour cohort100
            if strcmp(samples,'all'); obj.samples = size(EEG.data(1).series,1); %pour crises
            else; obj.samples = samples; end

            if any(strcmp(channels,ref)); channels(channels==ref)=[]; end
            idx = 1;
            %for ch = 1:size(EEG.data,2)                                %pour cohort100
            for ch = 1:size(EEG.data(1).series,2)                       %pour crises
                if any(strcmp(string(EEG.channels(ch).Name), channels))
                    obj.channels = [obj.channels, string(EEG.channels(ch).Name)];
                    %obj.rawData(:,idx) = EEG.data(start:start+samples-1,ch); %data normale
                    obj.rawData(:,idx) = EEG.data(1).series(start:start+samples-1,ch); 
                    % data crises
                    channels(channels==string(EEG.channels(ch).Name)) = [];
                    idx = idx +1;
                end
                if strcmp(string(EEG.channels(ch).Name), "A1")
                    %obj.mastoids = EEG.data(start:start+samples-1,ch:ch+1); %data normale
                    obj.mastoids = EEG.data(1).series(start:start+samples-1,ch:ch+1); 
                    % data crises
                end
            end

            obj.fs = EEG.channels.Fs;
        end

        function ppData = preprocess(obj,varargin)
            p = inputParser;
            addParameter(p, 'reref', '');
            addParameter(p, 'lowFilt', []);
            addParameter(p, 'highFilt', []);
            addParameter(p, 'window', [0,0]);
            addParameter(p, 'range', []);

            parse(p, varargin{:});
            rerefMode = p.Results.reref;
            lowFreqs = p.Results.lowFilt;
            highFreqs = p.Results.highFilt;
            windowSize = p.Results.window(1);
            overlap = p.Results.window(2);
            range = p.Results.range;
            
            ppData = obj.rereferencing(rerefMode);
            ppData = obj.filtering(ppData,lowFreqs,highFreqs);
            ppData = squeeze(obj.windowData(ppData,windowSize,overlap,range));
        end

        function rerefData = rereferencing(obj,mode)
            if strcmp(mode,"average")
                rerefData = obj.rereferenceToAverage;
            elseif strcmp(mode,"mastoids")
                rerefData = obj.rereferenceToMastoids;
            elseif strcmp(mode,"laplacian")
                rerefData = obj.rereferenceToLaplacian;
            else
                rerefData = obj.rawData;
            end
        end

        function rerefData = rereferenceToAverage(obj)
            ref = obj.rawData;
            if any(strcmp(obj.channels,"Fp1")); ref(:,obj.channels=="Fp1")=[]; end
            if any(strcmp(obj.channels,"Fp2")); ref(:,obj.channels=="Fp2")=[]; end
            if any(strcmp(obj.channels,"Cz")); ref(:,obj.channels=="Cz")=[]; end
            ref = mean(ref,2);
            rerefData = obj.rawData - ref;
        end

        function rerefData = rereferenceToMastoids(obj)
            ref = mean(obj.mastoids,2);
            rerefData = obj.rawData - ref;
        end

        function rerefData = rereferenceToLaplacian(obj)
            rerefData = zeros(size(obj.rawData));
            for i = 1:length(obj.laplace)
                surroundings = length(obj.laplace{i});
                if obj.laplace{i}(1) == "Fz" || obj.laplace{i}(1) == "Pz"
                    sum = (obj.rawData(:,obj.channels==obj.laplace{i}(2)) + obj.rawData(:,obj.channels==obj.laplace{i}(3)))/2;
                    for j = 4:surroundings
                        sum = sum + obj.rawData(:,obj.channels==obj.laplace{i}(j));
                    end
                    newRef = (1/(surroundings-2))*sum;
                else
                    sum = 0;
                    for j = 2:surroundings
                        if obj.laplace{i}(j) == "A1"; sum = sum + obj.mastoids(:,1);
                        elseif obj.laplace{i}(j) == "A2"; sum = sum + obj.mastoids(:,2);
                        else; sum = sum + obj.rawData(:,obj.channels==obj.laplace{i}(j));
                        end
                    end
                    newRef = (1/(surroundings-1))*sum;
                end
                rerefData(:,obj.channels==obj.laplace{i}(1)) = obj.rawData(:,obj.channels==obj.laplace{i}(1)) - newRef;
            end
        end

        function filtData = filtering(obj, data, lowFreqs, highFreqs)
            filtData = data(:,:,1);
            for band = 1:length(lowFreqs)
                [bHigh, aHigh] = butter(1,lowFreqs(band)/(obj.fs/2),"high"); 
                [bLow, aLow] = butter(2,highFreqs(band)/(obj.fs/2));
                filtHData = filtfilt(bHigh,aHigh,data);
                filtData(:,:,band) = filtfilt(bLow,aLow,filtHData);
            end
        end

        function windData = windowData(obj,data,window,overlap,range)
            if window == 0; windData = data;
            else
                window = window*obj.fs;
                overlap = overlap*obj.fs;
                windowStep = window - overlap;

                if ~isempty(range); data = data(range,:,:); end
                nbWindows = floor((size(data,1)-window)/windowStep)+1;
                windData = zeros([window size(data,2,3) nbWindows]);
                wIdx = 1;
                for w = 0:nbWindows-1
                    if w == 0; start = 1;
                    else; start = windowStep*w+1; end
                    wind = data(start:start+window-1,:,:);
                    if any(all(wind<0.0005 & wind>-0.0005,1),"all") || any(all(wind>500,1),"all") || any(all(wind<-500,1),"all")
                        windData(:,:,:,end) = [];
                    else
                        windData(:,:,:,wIdx) = wind;
                        wIdx = wIdx + 1;
                    end
                end 
            end
        end

    end
end
