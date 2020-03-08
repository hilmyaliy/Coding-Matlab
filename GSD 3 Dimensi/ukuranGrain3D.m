function [grainRadius,idxSingleGrain] = ukuranGrain3D(nGrain,imageGrainIdx,imSize,nTotalSize,grainNo,allGrainCentroid)

for iGrain = 1:nGrain

    idxSingleGrain= (imageGrainIdx(:) == grainNo(iGrain));
    imageGrain    = zeros(nTotalSize,1);
    imageGrain(idxSingleGrain) = 1;
    
        imageGrain = reshape(imageGrain,[imSize(1), imSize(2), imSize(3)]);

        [idxXX, idxYY, idxZZ] = find(imageGrainIdx == grainNo(iGrain));

        % principal component analysis to get principal direction (each column
        % = one principal component
        [coeff] = pca([idxXX idxYY idxZZ]);
        nVec = size(coeff,2);

        % Centroid
        X0 = allGrainCentroid(iGrain,1);
        Y0 = allGrainCentroid(iGrain,2);
        Z0 = allGrainCentroid(iGrain,3);

        for iVec = 1:nVec
            clear temp*

            % Obtain the value of secondary and tertiary direction
            a = coeff(1,iVec);
            b = coeff(2,iVec);
            c = coeff(3,iVec);

            [~,idxMax] = max([a,b,c]);

            if idxMax == 1
                tempX1 = [X0:imSize(1)]';
                tempX2 = [X0:-1:1]';
                nL1    = ones(length(tempX1),1);
                nL2    = ones(length(tempX2),1);
                tempY1 = nL1.*round(b./a.*(tempX1-X0) + Y0);
                tempY2 = nL2.*round(b./a.*(tempX2-X0) + Y0);    
                tempZ1 = nL1.*round(c./a.*(tempX1-X0) + Z0);
                tempZ2 = nL2.*round(c./a.*(tempX2-X0) + Z0);   
            elseif idxMax == 2
                tempY1 = [Y0:imSize(2)]';
                tempY2 = [Y0:-1:1]';
                nL1    = ones(length(tempY1),1);
                nL2    = ones(length(tempY2),1);
                tempX1 = nL1.*round(a./b.*(tempY1-Y0) + X0);
                tempX2 = nL2.*round(a./b.*(tempY2-Y0) + X0); 
                tempZ1 = nL1.*round(c./b.*(tempY1-Y0) + Z0);
                tempZ2 = nL2.*round(c./b.*(tempY2-Y0) + Z0);             
            elseif idxMax == 3
                tempZ1 = [Z0:imSize(3)]';
                tempZ2 = [Z0:-1:1]';
                nL1    = ones(length(tempZ1),1);
                nL2    = ones(length(tempZ2),1);
                tempX1 = nL1.*round(a./c.*(tempZ1-Z0) + X0);
                tempX2 = nL2.*round(a./c.*(tempZ2-Z0) + X0); 
                tempY1 = nL1.*round(c./c.*(tempZ1-Z0) + Y0);
                tempY2 = nL2.*round(c./c.*(tempZ2-Z0) + Y0);               
            end

            % Initialization of lines for measurement
            nLine1 = length(tempX1);
            nLine2 = length(tempX2);
            tempLine1 = zeros(nLine1,1);
            tempLine2 = zeros(nLine2,1);

            for iLine1 = 1:nLine1
                try % instead of checking boundary
                    tempLine1(iLine1,1) ...
                        = imageGrain(tempY1(iLine1),...
                                           tempX1(iLine1),...
                                           tempZ1(iLine1));
                end
            end

            for iLine2 = 1:nLine2
                try
                    tempLine2(iLine2,1) ...
                        = imageGrain(tempY2(iLine2),...
                                           tempX2(iLine2),...
                                           tempZ2(iLine2));
                end
            end

            % Find the boundary wihtin the lines
            bound1 = find(tempLine1 == 0, 1, 'first');
            if isempty(bound1)
                bound1 = nLine1;
            end

            bound2 = find(tempLine2 == 0, 1, 'first');
            if isempty(bound2)
                bound2 = nLine2;
            end

            % Calculate grain size to the boundary
            dX1 = tempX1(bound1) - tempX1(1);
            dY1 = tempY1(bound1) - tempY1(1);
            dZ1 = tempZ1(bound1) - tempZ1(1);
            dX2 = tempX2(bound2) - tempX2(1);
            dY2 = tempY2(bound2) - tempY2(1);
            dZ2 = tempZ2(bound2) - tempZ2(1);
            
            grainRadius(iGrain,2.*(iVec-1) + 1) ...
                = sqrt((dX1)^2 + (dY1)^2 + (dZ1)^2);
            grainRadius(iGrain,2.*(iVec-1) + 2) ...
                = sqrt((dX2)^2 + (dY2)^2 + (dZ2)^2);
                 
            grainAzimuth(iGrain,2.*(iVec-1) + 1) ...
                = acos(dZ1/grainRadius(iGrain,2.*(iVec-1) + 1));
            grainAzimuth(iGrain,2.*(iVec-1) + 2) ...
                = acos(dZ2/grainRadius(iGrain,2.*(iVec-1) + 2));
            
            grainInclination(iGrain,2.*(iVec-1) + 1) ...
                = atan(dY1/dX1)*180/pi;
            grainInclination(iGrain,2.*(iVec-1) + 2) ...
                = atan(dY2/dX2)*180/pi;

        end
end
idxNonZero          = all(grainRadius,2);
grainRadius         = grainRadius(idxNonZero,:);

 end



