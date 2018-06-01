function [smap,Color_cent,superlabel] = onescale_CB(i,image, seg_paras, bins, th, epsilon )
    % parameters used for multiple segmentations
    % bin of color histogram
    if ~exist('bins', 'var')
        bins = [8 16 16 4];       % CIE L*, a*, b*, and hue respectively
    end
    
    % threshold for merging adjacent superpixels
    if ~exist('th', 'var')
        th = 0.2;
    end
    
    % the parameter used in Eqn.(2)
    if ~exist('epsilon', 'var')
        epsilon = 0.1;
    end

 
    image = im2double(image);
    
    % compute quantize the image for color histogram generation
    Q = computeQuantMatrix(image, bins);
    
    smap = zeros(size(image,1), size(image,2));
    %addpath('./segment');
        % segment the image according to give parameters
        imsegs = im2superpixels(image, seg_paras);

        % for each region, compute the color histogram
        rh = computeRegionHist(Q, bins, imsegs.segimage);

        % to achieve better performance, merge adjacent regions if their
        % color distance is less then the given threshold
        imsegs2 = mergeAdjacentRegions_fast(rh, imsegs, th);
        % store superpixel .mat file(sigma=0.8) 
        superlabel=imsegs2.segimage;
        
        num_region = max(imsegs2.segimage(:));
        %fprintf('\t*** after merging, #num_region: %d\n', num_region);
       
        % compute the region color histogram for merged regions
        rh2 = computeRegionHist(Q, bins, imsegs2.segimage);

        % compute the color center for each region
        color_center = computeColorCenter(image, imsegs2.segimage);

        % for each pixel, compute its color distance to the color center of
        % the region which contains the pixel
        temp_color_weight = computeColorWeight(image, imsegs2.segimage, color_center, epsilon);

        % compute one superpixel-scale saliency map based on context
        % analysis
        temp_smap = computeOneScaleSmap_fast(i,rh2, imsegs2);
        
%         smaps{ix} = temp_smap;      
%         color_weight{ix} = temp_color_weight;        
        smap = (temp_smap - min(temp_smap(:))) / (max(temp_smap(:)) - min(temp_smap(:)) + eps);
        smap = uint8(smap * 255);
        Color_cent=temp_color_weight;
end
    
    
        