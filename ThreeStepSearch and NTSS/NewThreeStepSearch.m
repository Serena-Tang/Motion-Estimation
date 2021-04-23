% ECE 483 Digital Video Processing
% University of Victoria
% Jasmin Ford & Serena Tang
%
% This code was modified from the algorithm given by MATLAB
%
% Code begins below:
% Computes motion vectors using *NEW* Three Step Search method
%
% Based on the paper by R. Li, b. Zeng, and M. L. Liou
% IEEE Trans. on Circuits and Systems for Video Technology
% Volume 4, Number 4, August 1994 :  Pages 438:442
%
% Input
%   f2 : Target image frame
%   f1 : Reference (anchor) image frame
%   mbSize : Size of the macroblock
%   p : Search parameter  (read literature to find what this means)
%
% Ouput
%   u,v : motion vectors for predicting target frame frome reference frame
%   NTSScomputations: The average number of points searched for a macroblock
%
% Code modified by Jasmin Ford & Serena Tang from algorithm written by Aroh Barjatya

function [ u,v ] = NewThreeStepSearch( mbSize,p )

% f1 = imread('train01.tif'); %f1 = anchor frame
% f2 = imread('train02.tif'); %f2 = target frame

% f1 = imread('6.1.03.tiff'); %f1 = anchor frame
% f2 = imread('6.1.04.tiff'); %f2 = target frame

f1 = imread('motion05.512.tiff'); %f1 = anchor frame
f2 = imread('motion06.512.tiff'); %f2 = target frame


% display reference (anchor) image frame
figure
subplot(2,2,1)
imshow(f1)
hold on

% get dimensions of image frames
[height,width] = size(f1);
f1 = double(f1);
f2 = double(f2);

% create empty matrix to store motion vectors
motionVectors = zeros(2,height*width/mbSize^2); 

% create 3x3 matrix for keeping track of lowest MAD search point for a given
% search, initialize minimum value to 65537 (max value)
searchPoints = ones(3, 3) * 65537;

% we now take effectively log to the base 2 of p
% this will give us the number of steps required
L = floor(log10(p+1)/log10(2));   
stepMax = 2^(L-1);
computations = 0;

% we start off from the top left of the image
% we will walk in steps of mbSize
% for every marcoblock that we look at we will look for
% a close match p pixels on the left, right, top and bottom of it
mbCount = 1;
for i = 1 : mbSize : height-mbSize+1
    for j = 1 : mbSize : width-mbSize+1
        
        % the NEW three step search starts
        
        x = j;
        y = i;
        
        % In order to avoid calculating the center point of the search
        % again and again we always store the value for it from the
        % previous run. For the first iteration we store this value outside
        % the for loop, but for subsequent iterations we store the cost at
        % the point where we are going to shift our root.
        %
        % For the NTSS, we find the minimum first in the far away points
        % we then find the minimum for the close up points
        % we then compare the minimums and which ever is the lowest is where
        % we shift our root of search. If the minimum is the center of the
        % current window then we stop the search. If its one of the
        % immediate close to the center then we will do the second step
        % stop. And if its in the far away points, then we go doing about
        % the normal TSS approach
        % 
        % more details in the code below or read the paper/literature
        
        searchPoints(2,2) = calculateError(f1(i:i+mbSize-1,j:j+mbSize-1), ...
                                    f2(i:i+mbSize-1,j:j+mbSize-1),mbSize);
        stepSize = stepMax; 
        computations = computations + 1;
        % This is the calculation of the outer 8 points
        % m is row(vertical) index
        % n is col(horizontal) index
        % this means we are scanning in raster order
        for m = -stepSize : stepSize : stepSize        
            for n = -stepSize : stepSize : stepSize
                refBlkVer = y + m;   % row/Vert co-ordinate for ref block
                refBlkHor = x + n;   % col/Horizontal co-ordinate
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > height ...
                     || refBlkHor < 1 || refBlkHor+mbSize-1 > width)
                     continue;
                end
                costRow = m/stepSize + 2;
                costCol = n/stepSize + 2;
                if (costRow == 2 && costCol == 2)
                    continue
                end
                searchPoints(costRow, costCol ) = calculateError(f1(i:i+mbSize-1,j:j+mbSize-1), ...
                    f2(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                computations = computations + 1;
            end
        end
        
        % Now we find the vector where the cost is minimum
        % and store it ... 
        
        [dx, dy, min1] = determineMinError(searchPoints);      % finds which macroblock in imgI gave us min Cost
            
              
        % Find the exact co-ordinates of this point
        x1 = x + (dx-2)*stepSize;
        y1 = y + (dy-2)*stepSize;
            
        % Now find the costs at 8 points right next to the center point
        % (x,y) still points to the center
        
        stepSize = 1;
        for m = -stepSize : stepSize : stepSize        
            for n = -stepSize : stepSize : stepSize
                refBlkVer = y + m;   % row/Vert co-ordinate for ref block
                refBlkHor = x + n;   % col/Horizontal co-ordinate
                if ( refBlkVer < 1 || refBlkVer+mbSize-1 > height ...
                     || refBlkHor < 1 || refBlkHor+mbSize-1 > width)
                     continue;
                end
                costRow = m/stepSize + 2;
                costCol = n/stepSize + 2;
                if (costRow == 2 && costCol == 2)
                    continue
                end
                searchPoints(costRow, costCol ) = calculateError(f1(i:i+mbSize-1,j:j+mbSize-1), ...
                    f2(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                computations = computations + 1;
            end
        end
        
        % now find the minimum amongst this
        
        [dx, dy, min2] = determineMinError(searchPoints);      % finds which macroblock in imgI gave us min Cost
            
              
        % Find the exact co-ordinates of this point
        x2 = x + (dx-2)*stepSize;
        y2 = y + (dy-2)*stepSize;
        
        % the only place x1 == x2 and y1 == y2 will take place will be the
        % center of the search region
        
        if (x1 == x2 && y1 == y2)
            % then x and y still remain pointing to j and i;
            NTSSFlag = -1; % this flag will take us out of any more computations 
        elseif (min2 <= min1)
            x = x2;
            y = y2;
            NTSSFlag = 1; % this flag signifies we are going to go into NTSS mode
        else
            x = x1;
            y = y1;
            NTSSFlag = 0; % This value of flag says, we go into normal TSS
        end
        
        
        if (NTSSFlag == 1)
            % Now in order to make sure that we dont calcylate the same
            % points again which were in the initial center window we take
            % care as follows
            
            searchPoints = ones(3,3) * 65537;
            searchPoints(2,2) = min2;
            stepSize = 1;
            for m = -stepSize : stepSize : stepSize        
                for n = -stepSize : stepSize : stepSize
                    refBlkVer = y + m;   % row/Vert co-ordinate for ref block
                    refBlkHor = x + n;   % col/Horizontal co-ordinate
                    if ( refBlkVer < 1 || refBlkVer+mbSize-1 > height ...
                           || refBlkHor < 1 || refBlkHor+mbSize-1 > width)
                        continue;
                    end
                    
                    if ( (refBlkVer >= i - 1  && refBlkVer <= i + 1) ...
                            && (refBlkHor >= j - 1  && refBlkHor <= j + 1) )
                        continue;
                    end
                    
                    costRow = m/stepSize + 2;
                    costCol = n/stepSize + 2;
                    if (costRow == 2 && costCol == 2)
                        continue
                    end
                    searchPoints(costRow, costCol ) = calculateError(f1(i:i+mbSize-1,j:j+mbSize-1), ...
                         f2(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                    computations = computations + 1;
                end
            end
                
            % now find the minimum amongst this
        
            [dx, dy, min2] = determineMinError(searchPoints);      % finds which macroblock in imgI gave us min Cost
            
            % Find the exact co-ordinates of this point and stop
            x = x + (dx-2)*stepSize;
            y = y + (dy-2)*stepSize;            
            
        elseif (NTSSFlag == 0)
            % this is when we are going about doing normal TSS business
            searchPoints = ones(3,3) * 65537;
            searchPoints(2,2) = min1;
            stepSize = stepMax / 2;
            while(stepSize >= 1)  
                for m = -stepSize : stepSize : stepSize        
                    for n = -stepSize : stepSize : stepSize
                        refBlkVer = y + m;   % row/Vert co-ordinate for ref block
                        refBlkHor = x + n;   % col/Horizontal co-ordinate
                        if ( refBlkVer < 1 || refBlkVer+mbSize-1 > height ...
                            || refBlkHor < 1 || refBlkHor+mbSize-1 > width)
                            continue;
                        end
                        costRow = m/stepSize + 2;
                        costCol = n/stepSize + 2;
                        if (costRow == 2 && costCol == 2)
                            continue
                        end
                        searchPoints(costRow, costCol ) = calculateError(f1(i:i+mbSize-1,j:j+mbSize-1), ...
                                f2(refBlkVer:refBlkVer+mbSize-1, refBlkHor:refBlkHor+mbSize-1), mbSize);
                        computations = computations + 1;
                    
                    end
                end
        
                % Now we find the vector where the cost is minimum
                % and store it ... this is what will be passed back.
        
                [dx, dy, min] = determineMinError(searchPoints);      % finds which macroblock in imgI gave us min Cost
            
            
                % shift the root for search window to new minima point
                x = x + (dx-2)*stepSize;
                y = y + (dy-2)*stepSize;
            
                stepSize = stepSize / 2;
                searchPoints(2,2) = searchPoints(dy,dx);
            
            end
        end
        motionVectors(1,mbCount) = y - i;    % row co-ordinate for the vector
        motionVectors(2,mbCount) = x - j;    % col co-ordinate for the vector            
        
        current = 1;
        for k = 1:1:height/mbSize
            for l = 1:1:width/mbSize
                u(k,l) = motionVectors(1,current);
                v(k,l) = motionVectors(2,current);
                current = current + 1;
            end
        end
                
        mbCount = mbCount + 1;
        searchPoints = ones(3,3) * 65537;
    end
end

%NTSScomputations = computations/(mbCount - 1);

for i = 1:1:height/mbSize
    for j = 1:1:width/mbSize
        i1 = (i-1)*mbSize + 1;
        j1 = (j-1)*mbSize + 1;
        fp(i1:i1+mbSize-1,j1:j1+mbSize-1) = f2(u(i,j)+...
            i1:u(i,j)+i1+mbSize-1, v(i,j)+j1:v(i,j)+j1+mbSize-1);
    end
end

% construct error frame using predicted frame and reference frame
fe = fp - f1;

% calculate PSNR of error frame
psnr = 10*log10(255*255/mean(mean((fe).^2)));

% display anchor frame with motion vectors
[i,j] = meshgrid((mbSize+1)/2:mbSize:width,(mbSize+1)/2:mbSize:height);
quiver(i,j,u,v)     
title('Anchor frame with motion field')
hold off

% display motion field only
subplot(2,2,2)
u = -flipud(u);
v = flipud(v);
quiver(i,j,u,v)
axis([0 width 0 height]);
daspect([1 1 1]);
title(sprintf('Motion field with %d x %d search range',p,p))

% display predicted image frame with corresponding PSNR value
subplot(2,2,3)
imshow(uint8(fp))
title(sprintf('Predicted image (PSNR = %.4f)', psnr))

% display prediction error between predicted frame and anchor frame
subplot(2,2,4)
imshow(uint8(255 - abs(fe)))
title('Prediction error from New 3 Step Search')

end

% function to calculate error between current block being
% searched and reference block
function [ error ] = calculateError( currentBlock,refBlock,mbSize )

error = 0;

for i = 1:mbSize
    for j = 1:mbSize
        %when using MSE instead of MAD, psnr is improved
        error = error + abs((currentBlock(i,j) - refBlock(i,j)));
    end
end

%error = error / (mbSize^2);

end

% function to determine lowest error value of search points for a certain
% step
function [ dx,dy,minError ] = determineMinError( searchPoints )

[row,col] = size(searchPoints);

min = 65537;

for i = 1:row
    for j = 1:col
        if(searchPoints(i,j) < min)
            min = searchPoints(i,j);
            dx = j;
            dy = i;
        end
    end
end

minError = min;
end