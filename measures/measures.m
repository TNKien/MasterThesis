% Class Name:   measures
% Description:  Association measures
% Author:       Lise Cottin
% Date:         November 30, 2023

classdef measures
    properties
        nbIndicators
    end

    methods
        function obj = measures(nbIndicators)
            obj.nbIndicators = nbIndicators;
        end

        function matrix = association(obj,data,reref,fhbands,flbands,window,overlap)
            % Pre-processing
            ppData = data.preprocess('reref',reref,'window',[window,overlap]);

            % Connectivity matrix
            matrix = zeros([obj.nbIndicators length(fhbands) size(ppData,2,2)]);
            for i = 1:size(matrix,3)
                x = squeeze(ppData(:,i,:));
                parfor j = i+1:size(matrix,4)
                    y = squeeze(ppData(:,j,:));
                    matrix(:,:,i,j) = obj.measure(x,y,fhbands,flbands,data.fs);
                end
            end
        end

        function [xw,yw] = hannWindow(obj,x,y)
            x = (x - mean(x,1))./std(x,0,1);
            y = (y - mean(y,1))./std(y,0,1);
            hannWind = repmat(hann(size(x,1)),1,size(x,2));
            xw = hannWind.*x; 
            yw = hannWind.*y;
        end

    end
end
