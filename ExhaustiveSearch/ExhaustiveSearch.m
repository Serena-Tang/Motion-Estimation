% Working on modifying Pan's code from the assignment 4 solution, it works
% though!

% Input:  searchRange -- vector for search range for vert/horiz i.e. [8,8]
%         N -- uses N x N block size to search
% Output: u,v -- returns displacement vector
function [ u,v ] = ExhaustiveSearch( mbSize,searchRange )

% load reference and target image frames
% f1 = imread('train01.tif'); %f1 = anchor frame
% f2 = imread('train02.tif'); %f2 = target frame

% f1 = imread('6.1.03.tiff'); %f1 = anchor frame
% f2 = imread('6.1.04.tiff'); %f2 = target frame

f1 = imread('motion05.512.tiff'); %f1 = anchor frame
f2 = imread('motion06.512.tiff'); %f2 = target frame

% range_i = searchRange(1);
% range_j = searchRange(2);

% display reference (anchor) image frame
figure
subplot(2,2,1)
imshow(f1)
hold on

% get dimensions of image frames
[height,width] = size(f1);
f1 = double(f1);
f2 = double(f2);

computations = 0;
% initialize block size and search range
% N = 16;
% R = 32;
% MAD = 0;

% for each block in anchor frame (i = vertical, j = horizontal)
for i = 1:1:height/mbSize
    for j = 1:1:width/mbSize
        i1 = (i-1)*mbSize + 1;
        j1 = (j-1)*mbSize + 1;
        error = 1000000;
        % check vertical search range
        for k = max(1,i1-searchRange):1:min(height-mbSize+1,i1+searchRange)
            % check horizontal search range
            for l = max(1,j1-searchRange):1:min(width-mbSize+1,j1+searchRange)
                target_error = sum(sum(abs(f2(k:k+mbSize-1,l:l+mbSize-1)-...
                    f1(i1:i1+mbSize-1,j1:j1+mbSize-1))));
                if(target_error < error)
                    error = target_error;
                    u(i,j) = k - i1;
                    v(i,j) = l - j1;
                end
            end
        end
        % construct predicted image frame using determined motion vectors
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
title(sprintf('Motion field with %d x %d search range',searchRange,searchRange))

% display predicted image frame with corresponding PSNR value
subplot(2,2,3)
imshow(uint8(fp))
title(sprintf('Predicted image frame (PSNR = %.4f)', psnr))

% display prediction error between predicted frame and anchor frame
subplot(2,2,4)
imshow(uint8(255 - abs(fe)))
title('Prediction error from EBMA')
end