 classdef correlation < measures
    properties
        maxLag  {mustBeNumeric}
    end

    methods
        function obj = correlation(maxLag)
            obj@measures(2);
            obj.maxLag = maxLag;
        end

        function matrix = association(obj,data,reref,fhbands,flbands,window,overlap)
            % Pre-processing
            ppDataX = data.preprocess('reref',reref,'lowFilt',fhbands,'highFilt',flbands,'window',[window,overlap]);
            ppDataY = zeros([size(ppDataX,1)+(2*obj.maxLag) size(ppDataX,2,3) size(ppDataX,4)-2]);
            for w = 2:size(ppDataX,4)-1
                ppDataY(1:obj.maxLag,:,:,w-1) = ppDataX(end-obj.maxLag+1:end,:,:,w-1);
                ppDataY(obj.maxLag+1:obj.maxLag+size(ppDataX,1),:,:,w-1) = ppDataX(:,:,:,w);
                ppDataY(end-obj.maxLag+1:end,:,:,w-1) = ppDataX(1:obj.maxLag,:,:,w+1);
            end
%             disp(size(ppDataX)) ;
%             if size(ppDataX,4) > 2
                ppDataX(:,:,:,1) = []; ppDataX(:,:,:,end) = [];
%             end
%             disp(size(ppDataX)) ; disp(size(ppDataY)) ;
            % Connectivity matrix
            matrix = zeros([obj.nbIndicators size(ppDataX,3) size(ppDataX,2,2)]);
            for i = 1:size(matrix,3)
                x = squeeze(ppDataX(:,i,:,:));
                parfor j = i+1:size(matrix,4)
                    y = squeeze(ppDataY(:,j,:,:));
                    matrix(:,:,i,j) = obj.measure(x,y);                   
                end
            end
        end

        function corrInd = measure(obj,x,y)
            cc = obj.crossCorrelation(x,y);
            maxsCC = max(abs(cc),[],1);

            corr = obj.corrCrossCorrelation(cc);
            maxsCorr = max(abs(corr),[],1);

            corrInd = [mean(squeeze(maxsCC),2),mean(squeeze(maxsCorr),2)]';
        end

        % Cross-correlation
        function cc = crossCorrelation(obj,x,y)
            n = size(x,1);
            cc = zeros([2*obj.maxLag+1 size(x,2,3)]);
            meanx = mean(x,1);
            for lag = 1:2*obj.maxLag+1
                yc = y(lag:lag+n-1,:,:);
                meanyc = mean(yc,1);
                frac = 1./(std(x,0,1).*std(yc,0,1)*n);
                cross = zeros([1 size(x,2,3)]);
                for t=1:n
                    cross = cross + ((x(t,:,:)-meanx).*(yc(t,:,:)-meanyc));
                end
                cc(lag,:,:) = frac.*cross;
            end
        end

        % Corrected cross-correlation
        function corr = corrCrossCorrelation(obj,cc)
            corr = zeros([obj.maxLag size(cc,2,3)]);
            for lag = 1:obj.maxLag
                corr(lag,:,:) = (1/2).*(cc(lag+obj.maxLag+1,:,:)-cc(-lag+obj.maxLag+1,:,:));
            end
        end
    end
end
