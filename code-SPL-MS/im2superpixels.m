function imsegs = im2superpixels(im, seg_para)
    if nargin == 1
        seg_para = [0.8 100 100];
    end

%     prefix = num2str(floor(rand(1)*10000000));
%     fn1 = ['./tmpim' prefix '.ppm'];
%     fn2 = ['./tmpimsp' prefix '.ppm'];
%     segcmd = ['E:\playerkk\code\MATLAB\segment\segment ', num2str(seg_para(1)),... 
%         ' ', num2str(seg_para(2)), ' ', num2str(seg_para(3))];
% 
%     imwrite(im, fn1);
%     system([segcmd ' ' fn1 ' ' fn2]);
    if strcmp(class(im), 'uint8')
        im = double(im);
    end
    
    if max(im(:)) < 10
        im = double(im * 255);
    end
    
    segim = mexSegment(im, seg_para(1), seg_para(2), int32(seg_para(3)));
    imsegs = processSuperpixelImage(segim);
