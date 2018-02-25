%% clear
clc; clear all; close all;

%% path
addpath('./lib')
addpath('../vesselness3d/lib')

%% load image
im = imread3d('../vesselness3d/im/neuron.tif');

%% normalize
im = double(im); im = (im - min(im(:))) / (max(im(:)) - min(im(:)));

%% 3d vesselness
sigma = 3;
gamma = 2; 

imn = neuriteness3d(im,sigma,gamma);

%% plot
figure; imagesc(max(im,[],3)); colormap gray; 
set(gca,'ytick',[]); set(gca,'xtick',[]); axis image; axis tight;

figure; imagesc(max(imn,[],3)); colormap gray; 
set(gca,'ytick',[]); set(gca,'xtick',[]); axis image; axis tight;