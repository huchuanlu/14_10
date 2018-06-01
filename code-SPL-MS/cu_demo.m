clear all
clc
%matlabpool;
addpath('CBsegment');
sigma_g=1.5; % parameters for computing harris points
sigma_a=5;  % parameters for computing harris points
nPoints=30; % number of salient points 
thresh = 26; % for elimate the side point
inpath = '.\imgs\';
mkdir('out');
dir_im = dir([inpath '*.jpg']);  
tic;
for i =1:length(dir_im)
imName = dir_im(i).name; 
%imName='0_3_3317.jpg';
im=imread([inpath imName]);
w=size(im,1);
h=size(im,2);
%r = 400 / max(w, h); 
%im = imresize(im, r);
gray=im2double(rgb2gray(im));
%% detect and cut the meaningless frame of the image 
j=0;
for jj=1:16 % suppose the max size of the frame is 16
    box=gray(jj:w+1-jj,jj:h+1-jj);
    rim=[gray(jj,jj:h+1-jj),gray(w+1-jj,jj:h+1-jj),gray(jj:w+1-jj,jj)',gray(jj:w+1-jj,h+1-jj)'];
    if (max(rim)-min(rim))<0.3
       j=jj;
    end
end
if j>1
    I=zeros(w-2*j,h-2*j,3);
    for k=1:3
   image=im(:,:,k);
    I(:,:,k)=image(j+1:w-j,j+1:h-j);
    end
im=I;
end;
%%
im=im2double(im);
Mboost = BoostMatrix(im);
boost_im= BoostImage(im,Mboost);
[EnIm]= ColorHarris(boost_im,sigma_g,sigma_a,0.04);
[x_max,y_max,corner_im2,num_max]=getmaxpoints(EnIm,nPoints);
corner_im2 = elimatepoint(corner_im2,thresh); % elimate the points closing to the boundary of images
smap = multi_saliency(im,corner_im2 );
%% recover the original size of the saliency map
if(j>1) 
sm=zeros(w,h); 
sm(j+1:w-j,j+1:h-j)=smap;
else
sm=smap;
end
%sm = imresize(sm, [w h]);
%% save the map before filtering
sm=im2double(sm);
sm=(sm-min(sm(:)))/(max(sm(:))-min(sm(:))); 
imwrite(sm,['out\' imName]);
%% save the map after filtering
sm= guidedfilter(sm,sm,12,0.1);
sm=(sm-min(sm(:)))/(max(sm(:))-min(sm(:)));
sm=uint8(((sm-min(sm(:)))./(max(sm(:))-min(sm(:)))).*255);
imwrite(sm,['out\' imName]);
display(num2str(i));
%%
end
t=toc;