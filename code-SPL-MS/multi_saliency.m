function smap = multi_saliency(im,corner_im2 )
[w,h,~]=size(im);
t1 = zeros(w,h);
t2 = zeros(w,h);
a=[2,4,6,9,11]; 
for mm=1:length(a)
 k=a(mm);    
 [sal,Color_cent,PrH1,PrH0,out_PrH1,out_PrH0,ind]=cu_Saliency_map(k,im, corner_im2);
sal = im2double(sal);
sal = bayesian( sal, PrH1,PrH0,out_PrH1,out_PrH0,ind,w,h );
if(size(find(sal>0)) )
t1 = t1 + sal .* Color_cent;
t2 = t2 + Color_cent; 
end
end
smap = t1 ./ t2;
smap = uint8((smap - min(smap(:))) / (max(smap(:)) - min(smap(:)) + eps) * 255);
end

