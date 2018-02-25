function imn = neuriteness3d(im,sigma,gamma)
%%  neuriteness3d - vessel enhancement filtering
%   
%   REFERENCE:
%       C. Sazak and B. Obara,
%       Contrast-independent curvilinear structure enhancement in 3D
%       biomedical images, IEEE International Symposium on Biomedical 
%       Imaging. Melbourne, Australia, 1165-1168, 2017
%
%   INPUT:
%       im      - 3D gray image,
%       sigma   - Gaussian kernel sigma,
%       gamma   - parameter,
%
%   OUTPUT:
%       imn     - neuriteness
%
%   AUTHOR:
%       Boguslaw Obara

%% default parameters
if isempty(sigma);  sigma = 1;  end
if isempty(gamma);  gamma = 2;  end

%% normalize
im = double(im); im = (im - min(im(:))) / (max(im(:)) - min(im(:)));

%% convert image to grey-scale
% im = im2uint8(im); % I assume that Vesselness was defined for grey-scale images

%% convert image to double
% im = double(im);

%% second derivatives - Hessian
[hxx,hyy,hzz,hxy,hxz,hyz] = hessian3d(im,sigma);

%% normalized derivative - scale
hxx = power(sigma,gamma)*hxx;
hyy = power(sigma,gamma)*hyy;
hzz = power(sigma,gamma)*hzz;
hxy = power(sigma,gamma)*hxy;
hxz = power(sigma,gamma)*hxz;
hyz = power(sigma,gamma)*hyz;

%% eigen values and vectors
[l1,l2,l3] = ...
    eigen3d_m(hxx,hxy,hxz,hxy,hyy,hyz,hxz,hyz,hzz);

%% modified Hessian
alfa = -1/2;
l1p = l1 + alfa/2.*l2 + alfa/2.*l3;
l2p = l2 + alfa/2.*l1 + alfa/2.*l3;
l3p = l3 + alfa/2.*l1 + alfa/2.*l2;
l1 = l1p;
l2 = l2p;
l3 = l3p;

%% sort
[l1,~,~] = eigen_sort3d_m(l1,l2,l3);

%% neuriteness
lmin = min(l1(:));
imn = zeros(size(l1));
imn(l1<0) = l1(l1<0)./lmin;

end