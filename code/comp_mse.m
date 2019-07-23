function var = comp_mse(im)
im = im.^2;
[h,w] = size(im);
var = (1/(h*w-1))*sum(im(:));
end