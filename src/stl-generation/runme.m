% This script loads a stack of images of the vocal folds and converts them
% into an STL file

% Authors:        Emily Minch <evminch@wpi.edu>
%                 Rositsa Mihaleva <ramihaleva@wpi.edu>
%                 Alex Chiluisa <ajchiluisa@wpi.edu>
%                 Loris Fichera <lfichera@wpi.edu>
%
% Latest version: 6/14/2021

% Clean the workspace and close all the figures
clear, close all, clc

% What model are we working with?
modelID = 'L7'; % model ID
p = '\\research.wpi.edu\ROBOTICS\comet\NIDCD-2020\Vocal-Folds-Models-Bailly\LARYNX7\LARYNX7\L7a_pag_concaténé'; % path where the images are stored
prefix = 'L7a_pag_concat'; % image prefix
binarization_threshold = 29000;

% Define the path where the images are stored
addpath(p);

% Read the image size
i = imread(fullfile(p, [prefix '0000.tif']));
X = size(i,2);
Y = size(i,1);
Z = size(dir(fullfile(p, '*.tif')), 1);

% Define the region of interest
t = readtable('image-rois.csv', 'ReadRowNames', true);
roi = t({modelID},:);
xLB = roi.XLowerBound; xUB = roi.XUpperBound;
yLB = roi.YLowerBound; yUB = roi.YUpperBound;

% Define the sampling interval across the three axes
dZ = 2; dY = 2; dX = 2;

% Initialize the matrix where we are going to store all the images
VocalFolds = zeros((yUB - yLB)/dY, ...
                   (xUB - xLB)/dX,...
                   (Z/dZ));

% Load the images!
idz = 1;

f = waitbar(0,'Processing Images...');

for k = 0:dZ:(Z-1)
    idx = 1;
    idy = 1;
    
  % !FIXME the following code is not portable
  if k>9 && k<100
      myfilename = sprintf('L7a_pag_concat00%d.tif', k);
  elseif k>99 && k<1000
       myfilename = sprintf('L7a_pag_concat0%d.tif', k);
  elseif k>999 && k<10000
       myfilename = sprintf('L7a_pag_concat%d.tif', k);
  else
      myfilename = sprintf('L7a_pag_concat000%d.tif', k);
  end

  raw_image = imread(myfilename);             % Read the image
  cropped_image = raw_image(yLB:yUB,xLB:xUB); % Crop the image
  downscaled_image = zeros((yUB-yLB)/dY,(xUB-xLB)/dX); % Initialize the downscaled image
 
  for p = 1 : dX : (xUB-xLB)
      for q = 1 : dY : (yUB-yLB)
          downscaled_image(idx,idy) = cropped_image(q,p);
          idx=idx+1; 
      end
      idx = 1;
      idy = idy+1;
  end
  
  % Binarize the image
  [L, W] = size(downscaled_image);
  BW = zeros(L,W);  
  BW = downscaled_image > binarization_threshold;
              
  VocalFolds(:,:,idz) = BW;                                       % match dimensions
  idz = idz + 1;
  waitbar(k/(Z-1),f,'Processing Images...');
end

close(f)

make_STL_of_Array('L7aV3.stl',VocalFolds,...
                  VoxelSizeX * dX,...
                  VoxelSizeY * dY,...
                  VoxelSizeZ * dZ);