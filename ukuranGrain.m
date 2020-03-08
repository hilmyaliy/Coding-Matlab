function [grainRadius,idxSingleGrain,coeff,idxXX, idxYY] = ukuranGrain(nGrain,grainNo,nTotalSize,imsize,allGrainCentroid,water)

for iGrain = 1:nGrain
   
    idxSingleGrain      = (water(:) == grainNo(iGrain));
    imageGrain    = zeros(nTotalSize,1);
    imageGrain(idxSingleGrain) = 1;
    
    
    %if imsize(3) == 1 % 2D-------------------------------------------------
        imageGrain    = reshape(imageGrain,[imsize(2), imsize(1)]);
        [idxXX, idxYY] = find(water == grainNo(iGrain)); %menemukan posisi grain
        
        % Obtain principal direction
        [coeff] = pca([idxXX idxYY]);
        
        nVec = size(coeff,2);
        X0 = allGrainCentroid(iGrain,1); %pusat
        Y0 = allGrainCentroid(iGrain,2); %pusat

        for iVec = 1:nVec
            clear temp*

            % Find index of a straight line in both direction on principal axis
            a = coeff(1,iVec);
            b = coeff(2,iVec);

            if a >= b %m=dy/dx
                tempX1 = [X0:imsize(2)]';
                tempX2 = [X0:-1:1]';
                
                nL1    = ones(length(tempX1),1);
                nL2    = ones(length(tempX2),1);
                tempY1 = nL1.*round(b./a.*(tempX1-X0) + Y0); %b/a adalah slope
                tempY2 = nL2.*round(b./a.*(tempX2-X0) + Y0);    
                
            elseif b > a %m=dx/dy
                tempY1 = [Y0:imsize(1)]';
                tempY2 = [Y0:-1:1]';
                nL1    = ones(length(tempY1),1);
                nL2    = ones(length(tempY2),1);
                tempX1 = nL1.*round(a./b.*(tempY1-Y0) + X0);
                tempX2 = nL2.*round(a./b.*(tempY2-Y0) + X0);        
            end

            % Initialization of lines for measurement
            nLine1 = length(tempX1);
            nLine2 = length(tempX2);
            tempLine1 = zeros(nLine1,1);
            tempLine2 = zeros(nLine2,1);

            for iLine1 = 1:nLine1
                try % instead of checking boundary
                tempLine1(iLine1,1) = imageGrain(tempY1(iLine1),...
                                                       tempX1(iLine1));
                end
            end

            for iLine2 = 1:nLine2
                try
                tempLine2(iLine2,1) = imageGrain(tempY2(iLine2),...
                                                       tempX2(iLine2));
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
            
            dX2 = tempX2(bound2) - tempX2(1);
            dY2 = tempY2(bound2) - tempY2(1);
            
            % Compute 2 radius from PCA at a time
            grainRadius(iGrain,2.*(iVec-1) + 1)= sqrt(dX1^2 + dY1^2); %1 dan 2
            grainRadius(iGrain,2.*(iVec-1) + 2)= sqrt(dX2^2 + dY2^2); %3 dan 4

            grainAzimuth(iGrain,2.*(iVec-1) + 1) = atan(dX1/dY1)*180/pi;
            grainAzimuth(iGrain,2.*(iVec-1) + 2) = atan(dX2/dY2)*180/pi;
            
                %menampilkan grafik diameter
                %if imsize(3) == 1 %jika 2 dimensi
                hold on
                plot(tempX1(1:bound1),tempY1(1:bound1),'r','LineWidth',2)
                plot(tempX2(1:bound2),tempY2(1:bound2),'r','LineWidth',2)
                %end           
        end
    %end
end


idxNonZero          = all(grainRadius,2);

grainRadius         = grainRadius(idxNonZero,:);



end



