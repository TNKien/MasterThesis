classdef phase < measures
    methods
        function obj = phase
            obj@measures(3);
        end

        function phaseInd = measure(obj,x, y,fhbands,flbands,fs)
            [plv, pli, wpli, f] = obj.phaseSync(x,y,fs);

            phaseInd = zeros(obj.nbIndicators, length(fhbands));
            for b = 1:length(fhbands)
                phaseInd(1,b) = mean(plv(f>=fhbands(b) & f<flbands(b),:),1);
                phaseInd(2,b) = mean(pli(f>=fhbands(b) & f<flbands(b),:),1);
                phaseInd(3,b) = mean(wpli(f>=fhbands(b) & f<flbands(b),:),1);
            end
        end

        % PLV, PLI and wPLI
        function [plv,pli,wpli,f] = phaseSync(obj,x,y,fs)
            [xw,yw] = obj.hannWindow(x,y);
            [sxy,f] = cpsd(xw,yw,[],[],[],fs);

            plv = abs(mean(exp(1i.*angle(sxy)),2));
            pli = abs(mean(sign(imag(sxy)),2));
            num = abs(mean(imag(sxy),2));
            den = mean(abs(imag(sxy)),2);
            wpli = num./den;
        end

    end
end
