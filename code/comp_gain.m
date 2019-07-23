function [G, sub] = comp_gain(dec_im, lev, type)
%computes coding gain: assumes size of the image to be a power of 2

a_mean = 0;
sub_var = 1;
[h, w] = size(dec_im);

if strcmp(type,'uniform')
    %extract each subband
    hor = ones(1,2^lev)*w/2^lev;
    vert = ones(1,2^lev)*h/2^lev;
    sub = mat2cell(dec_im, vert, hor);
    M = 4^lev;
    %compute geometric mean
    for i = 1:M
        a_mean = a_mean + var(double(sub{i}(:)));
        sub_var = sub_var*(var(double(sub{i}(:))).^(1/M));
    end
    a_mean = a_mean/M;
elseif strcmp(type,'pyramidal')
    l = lev;
    sub = [];
    M = 3*lev+1;
    %extract subbands
    while ~isequal(l,0)
        [~, aux{l}] = comp_gain(dec_im(1:h/(2^(lev-l)),1:w/(2^(lev-l))),1, 'uniform');
        if l ~= 1  %discard the upper left square and mantain the others
            sub = cat(2,sub,aux{l}(2:4));
        else  % at the last iteration keep all the squares
            sub = cat(2,sub,aux{l}(1:4));
        end
        l = l-1;
    end
    %compute geometric mean
    for i = 1:M
        a_mean = a_mean + var(double(sub{i}(:)))*(numel(sub{i})/(w*h));
        sub_var = sub_var*(var(double(sub{i}(:))).^(numel(sub{i})/(w*h)));
    end
else
    error('unsupported decompostion strategy');
end

%compute gain
%G = var(double(im(:)))/sub_var;
G = a_mean/sub_var;
end