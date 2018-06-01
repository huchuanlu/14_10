function saliencyMap = computeOneScaleSmap_fast(i,regionHist, imsegs)
    segimage = imsegs.segimage;
    [h w c] = size(segimage);
    
    [x y] = meshgrid(-w/2+1:w/2, -h/2+1:h/2); 
    x = x.^2;
    y = y.^2;
    
    adjmat = imsegs.adjmat;
    
    num_region = max(imsegs.segimage(:));
    
    saliencyMap = zeros(h, w);
    
    edge_num=zeros(num_region,1);
    
    spstats = regionprops(segimage, 'all');
    area = zeros(num_region, 1);
    spatial_prior = zeros(num_region, 1);
    c = 3;
    for ix = 1 : num_region
        area(ix) = length(spstats(ix).PixelIdxList);
        adjmat(ix, ix) = 0;
        
        indxy=spstats(ix).PixelList;
        edge_num(ix)=length(find(indxy(:,1)==3))+length(find(indxy(:,1)==(size(image,1)-2)))+length(find(indxy(:,2)==3))+length(find(indxy(:,2)==(size(image,2)-2)));
        
        temp_x_dist = mean(x(spstats(ix).PixelIdxList));
        temp_y_dist = mean(y(spstats(ix).PixelIdxList));
        spatial_prior(ix) = exp(-c^2*temp_x_dist/w^2 - c^2*temp_y_dist/h^2);
    end
    
    area_weight = repmat(area', [num_region, 1]) .* double(adjmat);
    
    area_weight = area_weight ./ repmat(sum(area_weight, 2), [1, num_region]);
    
    node_weight = area / h /w / 0.52;
    node_weight = 1 ./ (1 + (node_weight).^11);
    
    region_dist = zeros(num_region, num_region);
    
    ind = find(adjmat);
    for ix = 1 : length(ind)
        [x y] = ind2sub([num_region, num_region], ind(ix));
        region_dist(x, y) = histDist(regionHist(x,:), regionHist(y,:));
    end
   % spatial_prior=1;
    temp_saliency = -log(1 - region_dist) .* area_weight;
    saliency = spatial_prior .* sum(temp_saliency, 2) .* node_weight;
   %saliency=spatial_prior;
    th1=70;
    th2=100;
    if(i<5)
        th1=50;
        th2=80;
    end
    % th1=100;
    %th2=200;
    %if(i<5)
    %    th1=80;
    %    th2=150;
    %end
    
    saliency=saliency.*exp(-edge_num/th1);
    for ix = 1 : num_region
        if(edge_num(ix)>th2)
            saliency(ix)=0;
       end
    end
    
    for ix = 1 : num_region
        saliencyMap(spstats(ix).PixelIdxList) = repmat(saliency(ix), size(spstats(ix).PixelIdxList));
    end
    
    %saliencyMap = saliencyMap .* spatialPrior;