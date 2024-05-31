classdef information < measures
    properties
        bins
        % tau     = [32,16,8,4,16,32];
        tau     = [128,64,32,16,32,16]./1000;   % seconds
        kernel  = 3;
    end

    methods
        function obj = information(bins)
            obj@measures(2);
            obj.bins = bins;
        end

        function matrix = association(obj,data,reref,fhbands,flbands,window,overlap)
            % Pre-processing
            ppData = data.preprocess('reref',reref,'lowFilt',fhbands,'highFilt',flbands,'window',[window,overlap]);

            % Connectivity matrix
            matrix = zeros([obj.nbIndicators size(ppData,3) size(ppData,2,2)]);
            for i = 1:size(matrix,3)
                x = squeeze(ppData(:,i,:,:));
                parfor j = i+1:size(matrix,4)
                    y = squeeze(ppData(:,j,:,:));
                    matrix(:,:,i,j) = obj.measure(x,y,data.fs);
                end
            end
        end

        function infInd = measure(obj,x,y,fs)
            x = (x - mean(x,1))./std(x,0,1);
            y = (y - mean(y,1))./std(y,0,1);
            infInd = zeros(obj.nbIndicators,size(x,2));
            for b = 1:size(x,2)
                mi = zeros(1,size(x,3)); 
                wsmi = zeros(1,size(x,3));
                for w = 1:size(x,3)
                    xw = squeeze(x(:,b,w)); 
                    yw = squeeze(y(:,b,w));
                    mi(w) = obj.mutualInformation(xw,yw,[],[]);
                    wsmi(w) = obj.wSMI(xw,yw,obj.tau(b),fs);
                end
                infInd(:,b) = [mean(mi,2),mean(wsmi,2)];
            end
        end

        % MI
        function mi = mutualInformation(obj,x,y,weights,bins)
            if isempty(weights); weights = ones(obj.bins,obj.bins); end
            if isempty(bins); bins = obj.bins; end

            edges = linspace(min([min(x),min(y)]), max([max(abs(x)),max(abs(y))+0.1]), bins+1);

            hxy = histcounts2(x, y, edges, edges); pxy = hxy/sum(hxy,"all");
            hx = histcounts(x, edges); px = hx/sum(hx);
            hy = histcounts(y, edges); py = hy/sum(hy);
            
            mi = 0;
            for i = 1:bins
                for j = 1:bins
                    if pxy(i,j) > 0
                        mi = mi + weights(i,j)*pxy(i,j)*log10(pxy(i,j)/(px(i)*py(j)));
                    end
                end
            end
        end

        % wSMI
        function wsmi = wSMI(obj,x,y,t,fs)
            tauPts = t*fs;
            % Symbolic transform
            xHat = obj.symbolicTransform(x,tauPts);
            yHat = obj.symbolicTransform(y,tauPts);
            
            % weighted Mutual information
            weights = [0,1,1,1,0,1;1,0,1,1,1,0;1,1,0,0,1,1;1,1,0,0,1,1;0,1,1,1,0,1;1,0,1,1,1,0];
            wsmi = obj.mutualInformation(xHat,yHat,weights,factorial(obj.kernel));
            wsmi = (1/log10(factorial(obj.kernel))) * wsmi;
        end

        function symbolic = symbolicTransform(obj,signal,tau)
            symbolic = zeros(length(signal)-tau*(obj.kernel-1),1);
            for s = 1:length(symbolic)
                range = s:tau:s+(obj.kernel-1)*tau;
                if (signal(range(2))<signal(range(1))) && (signal(range(2))>signal(range(3)))
                    symbolic(s) = 1;    % negative line
                elseif (signal(range(3))<signal(range(1))) && (signal(range(3))>signal(range(2)))
                    symbolic(s) = 2;    % U left
                elseif (signal(range(1))<signal(range(2))) && (signal(range(1))>signal(range(3)))
                    symbolic(s) = 3;    % bridge right
                elseif (signal(range(1))<signal(range(3))) && (signal(range(1))>signal(range(2)))
                    symbolic(s) = 4;    % U right
                elseif (signal(range(2))<signal(range(3))) && (signal(range(2))>signal(range(1)))
                    symbolic(s) = 5;    % positive line
                elseif (signal(range(3))<signal(range(2))) && (signal(range(3))>signal(range(1)))
                    symbolic(s) = 6;    % bridge left
                end
            end
        end

    end
end
