% Create an image stack of vocal fold images
% Turn them into a 3D file

% NOTE: Takes a long ass time, do not do anything else during calculations
% Need to use matlab server, 16G ram don't cut it

addpath('\Vocal-Folds-Models\LARYNX1\LARYNX1\SlicesZ');
%  idx = 0:141;
%  # image file paths
%  ipaths =cellfun(@(x) sprintf('IMG_ %04 d . JPG',x),idx,'uni',false);
%  # read images into cell array
%  I = cellfun(@imread, ipaths,'uni',false);
%  # Convert to gray-scale 
%  Igray = cellfun(@rgb2gray,I,'uni',false);
%  # Concatenate 3 gray-scale images into single 3D matrix 
%  myImage = cat(3,Igray{:});
clear;clc;

numfiles = 1301;
% mydata = cell(1, numfiles);
dZ =5; %1 per how many slices in the z direction
dx = 2; %1 per how many slices in the x
dy = 2; %1 per how many slices in the y
LB = 570;
RB = 1760;
length = RB-LB;
VocalFolds = zeros(1340/dy,length/dx,(1760/dZ));
count = 1;

for k = 200:dZ:(numfiles-1)
    countx = 1; 
    county = 1;
    
  if k>9 && k<100
      myfilename = sprintf('000% d . tif', k);
  elseif k>99 && k<1000
       myfilename = sprintf('00% d . tif', k);
  elseif k>999 && k<10000
       myfilename = sprintf('0% d . tif', k);
  elseif k>9999
       myfilename = sprintf('% d . tif', k);
  else
      myfilename = sprintf('0000% d . tif', k);
  end
      
  % myfilename = sprintf('file % d . txt', k);
  % mydata{k} = importdata(myfilename);
  imagek_d = imread(myfilename);
  imagek_n = imagek_d(1:1340,LB+1:RB);
  imagek = zeros(1340/dy,length/dx);
  for p = 1:dx:length 
      for q = 1:dy:1340 
          imagek(countx,county) = imagek_n(q,p);
          countx=countx+1; 
      end
      countx=1;
      county = county+1;
  end
  % BW = im2bw(imagek,0.00001); % stl writer only works with BW image
  [L W] = size(imagek);
  BW = zeros(L,W);
  for i = 1:L
      for j = 1:W
          if imagek(i,j) > 1284
              BW(i,j) = 1;
          end
      end
  end
              
  % imagek=rgb2gray(imagek);
  % VocalFolds(:,:,count) = imagek; 
  VocalFolds(:,:,count) = BW;
  count=count+1;
  end

% run make_stl after, it saves RAM an dlets you do larger matrices
% make_STL _of _Array('L1_v2 . stl',VocalFolds,0.045142*dx,0.045142*dy,(0.045142*dZ));

% save('Larynx1','VocalFolds')  % File too large