function [sal_super,Color_cent,PrH1,PrH0,out_PrH1,out_PrH0,ind]=cu_Saliency_map(i,im,corner_im)
input_im = rgb2lab(im);%LAB color space
      seg_paras = [ 0   100 1200;
                    0.2 150 1200;
                    0.3 175 1200;
                    0.4	200	1200;
                    0.5 225 1200;
                    0.6	250	1200;
                    0.7 275 1200;
                    0.8	300	1200;
                    0.9 325 1200;
                    1.0 350 1200;
                    1.2 400 1200;
                    1.4 450 1200];
  [row,col] = size(corner_im);
  [y,x] = ind2sub([row,col],find(corner_im == 1));
  if isempty(y) || length(y)<=2
     corner_im(floor(row/4):floor(row*3/4),floor(col/4):floor(col*3/4))=1; 
     [y,x] = ind2sub([row,col],find(corner_im == 1));
  end
  dt = DelaunayTri(x,y);
  if(~size(dt,1))
      return;
  end
  [k,av] = convexHull(dt);
[sal_super,Color_cent,sulabel_im] = onescale_CB(i,im, seg_paras(i,:));
sal_super=im2double(sal_super);
 BW = roipoly(corner_im,x(k),y(k));
 pixel = regionprops(BW,'all');
 ind = pixel.PixelIdxList;
 out_ind = setdiff(1:row*col,ind);
 ind0=ind;
%% 
STATS = regionprops(sulabel_im, 'all');
sup_num = numel(STATS);
for r = 1:sup_num %% check the superpixel not along the image sides
    indxy = STATS(r).PixelIdxList;
  if (numel(intersect(indxy,ind)) > 0.60 * numel(indxy))
        ind = unique([ind;indxy]);
        out_ind = setdiff(out_ind,indxy');
    else
       out_ind = unique([out_ind,indxy']);
       ind = setdiff(ind,indxy);
    end
end
avbi=length(ind)/av;
if (  avbi<0.25 ||avbi>1.4 ||isempty(ind) || isempty(out_ind)||(mean(sal_super(ind0))-mean(sal_super(ind)))>0.15)
    ind=ind0;
end
%%
corner=zeros(row,col);
corner(ind)=1;
thresh=15;
corner = elimatepoint(corner,thresh);
corner=bwlabel(corner);
num_re=max(corner(:));
if (num_re>1&&num_re<4)
    id1=find(corner==1);
    for k=2:num_re
        idk=find(corner==k);
        if(length(idk)/length(id1)>1.6 )
            corner(id1)=0;
            id1=idk;
        else if(length(id1)/length(idk)>1.6 )
                corner(idk)=0;
            else if(mean(sal_super(idk))-mean(sal_super(id1)))>0.3
                    corner(id1)=0;
                   id1=idk;
                else if(mean(sal_super(id1))-mean(sal_super(idk)))>0.3
                     corner(idk)=0;   
                    end
                end
            end
        end
    end
end  
ind=find(corner);
%%
 if ( length(ind)/av<0.25 )
    ind=ind0;
 end
 out_ind=setdiff(1:row*col,ind);

mat_im = [reshape(input_im(:,:,1),1,row*col);reshape(input_im(:,:,1),1,row*col);reshape(input_im(:,:,1),1,row*col)];
maxValO = max(mat_im(:,ind),[],2);
minVal0 = min(mat_im(:,ind),[],2);
maxValB = max(mat_im(:,out_ind),[],2);
minValB = min(mat_im(:,out_ind),[],2);
numBin=[60,60,60]; % Number of bins in histogram (if 2D histogram are used, numBin=[numBin1, numBin2];)
smoothFactor=[5,6,6]; % Smoothing factor
smoothingKernel=cell(1,3);
% PrH1 = 0.2;%numel(ind)/(row*col)
% PrH0 = 1-PrH1;
% out_PrH1 =0.2;
% out_PrH0 = 1-out_PrH1;
PrH1 = 1;%numel(ind)/(row*col)
PrH0 = 1;
out_PrH1 =1;
out_PrH0 = 1;

for i = 1:3
    cur_im = input_im(:,:,i);
    dataMat = cur_im(ind);
    [innerHist,innerBin] = ComputeHistogram_(dataMat,numBin(i),minVal0(i),maxValO(i));
    smoothingKernel{i}=getSmoothKernel_(smoothFactor(i));
    innerHist=filterDistribution_(smoothingKernel{i},innerHist',numBin(i));
    
    dataMat = cur_im(out_ind);
    [outerHist,outerBin] = ComputeHistogram_(dataMat,numBin(i),minValB(i),maxValB(i));
    smoothingKernel{i}=getSmoothKernel_(smoothFactor(i));
    outerHist=filterDistribution_(smoothingKernel{i},outerHist',numBin(i));
    
    PrO_H1 = innerHist(innerBin);
    PrO_H0 = outerHist(innerBin);
    PrH1=PrH1.*PrO_H1;
    PrH0=PrH0.*PrO_H0;
    
    
    PrB_H1 = innerHist(outerBin);
    PrB_H0 = outerHist(outerBin);
    out_PrH1=out_PrH1.*PrB_H1;
    out_PrH0=out_PrH0.*PrB_H0;
end


function [intHist,binInd] = ComputerColorHist(dataMat,color_bin)
%dataMat: innerIm: inner ind image
%         outerIm: outer ind image
%intHist: the output histgram for the input im 
%binInd : bin of each number of the input im
% size = size(dataMat,2);
numBin = color_bin^3;
dataMat_max = max(dataMat(:));
dataMat_min = min(dataMat(:));
dataMat = (dataMat - dataMat_min)/(dataMat_max - dataMat_min);% 归一化
binInd = floor(dataMat(1,:) * color_bin * color_bin + dataMat(2,:) * color_bin + dataMat(3,:)) + 1;% ind1-4096
intHist=zeros(numBin,1);% 
for i = 1:length(dataMat)
    intHist(binInd(i))=intHist(binInd(i))+1;
end


function [intHist,binInd]=ComputeHistogram_(dataMat,numBin,minVal,maxVal)
% currently only 1 histograms
% dataMat:L orA orB with inner or outer index

%output:
%  intHist: 
%  binInd :

%binInd=ones(length(dataMat));
binInd=max( min(ceil(numBin*(double(dataMat-minVal)/(maxVal-minVal))),numBin),1);% 将像素值 转换到[1 - 60]之间，大小中原图等同
intHist=zeros(numBin,1);
for i = 1:length(dataMat)
    intHist(binInd(i))=intHist(binInd(i))+1;
end
%% Get smoothing filter
function [smKer]=getSmoothKernel_(sigma)

if sigma==0
    smKer=1;
    return;
end

dim=length(sigma); % row, column, third dimension
sz=max(ceil(sigma*2),1);
sigma=2*sigma.^2;

if dim==1
    d1=-sz(1):sz(1);
    
    smKer=exp(-((d1.^2)/sigma));
    
elseif dim==2
    [d2,d1]=meshgrid(-sz(2):sz(2),-sz(1):sz(1));
    
    smKer=exp(-((d1.^2)/sigma(1)+(d2.^2)/sigma(2)));
    
elseif dim==3
    [d2,d1,d3]=meshgrid(-sz(2):sz(2),-sz(1):sz(1),-sz(3):sz(3));
    
    smKer=exp(-((d1.^2)/sigma(1)+(d2.^2)/sigma(2)+(d3.^2)/sigma(3)));
    
else
    error('Not implemented');
end

smKer=smKer/sum(smKer(:));





%% Smooth distribution
function dist=filterDistribution_(filterKernel,dist,numBin)

if numel(filterKernel)==1
    dist=dist(:)/sum(dist(:));
    return;
end

numDim=length(numBin);

if numDim==1
    %smoothDist=conv(dist,filterKernel,'same');
    
    lenDist=length(dist);
    hlenKernel=(length(filterKernel)-1)/2;

    dist=[dist(1)*ones(1,hlenKernel),dist,dist(end)*ones(1,hlenKernel)];
    dist=conv(dist,filterKernel);
    lenSmoothDist=length(dist);
    offset=(lenSmoothDist-lenDist)/2;
    dist=dist((offset+1):(lenSmoothDist-offset));
    
elseif numDim==2
    dist=reshape(dist,numBin);
    
    dist=conv2(filterKernel,filterKernel,dist,'same');
    
else
    dist=reshape(dist,numBin);
    
    for i=1:numDim
        fker=ones(1,numDim);
        fker(i)=length(filterKernel);
        fker=zeros(fker);
        fker(:)=filterKernel(:);
        
        dist=convn(dist,fker,'same');
        
    end
    
end

dist=dist(:)/sum(dist(:));
    





 