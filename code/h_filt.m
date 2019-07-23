function [im_lp, im_hp] = h_filt(im)
%computes the decomposition row-wise


[h, w] = size(im);

%compute wavelet filters coefficients
[Lp,Hp,~,~] = wfilters('haar');

%initialize quantities
L = numel(Lp);
im_lp = zeros(h,w+3*L-1);
im_hp = im_lp;

%pad the image along rows and filter
for i = 1:h
    im_lp(i,:) = conv(double(padarray(im(i,:)', length(Lp),'circular'))',Lp);
    im_hp(i,:) = conv(double(padarray(im(i,:)', length(Hp),'circular'))',Hp);
end

%take the correct values and downsample along rows
im_lp = rot90(downsample(rot90(im_lp(:,2*L:w + 2*L-1),3),2));
im_hp = rot90(downsample(rot90(im_hp(:,2*L:w + 2*L-1),3),2));
end