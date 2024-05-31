classdef coherence < measures
    methods
        function obj = coherence
            obj@measures(2);
        end

        function cohInd = measure(obj,x,y,fhbands,flbands,fs)
            [cohy, f] = obj.coherency(x,y,fs);

            cohInd = zeros(obj.nbIndicators, length(fhbands));
            for b = 1:length(fhbands)
                meanCohy = mean(cohy(f>=fhbands(b) & f<flbands(b),:),1);
                cohInd(:,b) = [mean(abs(meanCohy),2), mean(imag(meanCohy),2)];
            end
            
        end

        % Coherency
        function [coherency, f] = coherency(obj,x,y,fs)
            [xw,yw] = obj.hannWindow(x,y);
            [sxy,f] = cpsd(xw,yw,[],[],[],fs);
            [sxx,~] = cpsd(xw,xw,[],[],[],fs);
            [syy,~] = cpsd(yw,yw,[],[],[],fs);
            coherency = sxy./sqrt(sxx.*syy);
        end
    end
end
