function [quant, encoded, dictionary, len] = quant_subband(dimg, q_table)
%works only for uniform decomposition, computes the quantization and
%the following entropy coding (huffman):
% quant: is the decomposed image quantized
% encoded: is the decomposed image quantized and entropy coded
% dictionary: is the hufffman dictionary for each subband
% len: is the length of huff. dict for each subband

[h, w] = size(dimg);
lev = log2(length(q_table))/2;      %it should be a log4
M = 4^lev;                         %number of subbands
h_step = w/2^(lev);               %width in pixel of a subband
v_step = h/2^(lev);                 %height in pixel of a subband
quant = zeros(h,w);
len = zeros(1,M);
encoded = []; dictionary = {};

%extracts subbands
hor = ones(1,2^lev)*w/2^lev;
vert = ones(1,2^lev)*h/2^lev;
subb = mat2cell(dimg, vert, hor);

%quantize each subband
for i = 1:M
    if q_table(i) == 0  %if no bit for that subband discard the values and don't code anything
        quant((v_step*(mod(i-1,2^lev))+1:v_step*(mod(i-1,2^lev))+v_step),...
            (h_step*ceil(i/2^lev-1)+1:h_step*ceil(i/2^lev-1)+h_step)) = zeros(v_step,h_step);
        dictionary{i} = 0;
        len(i) = len(i-1);  %ok because in practice I can't go here at
        %the first iteration otherwise I'm not coding anything
    else
        a = reshape(subb{i},[1,numel(subb{i})]);
        delta = (max(a)-min(a))/2^(q_table(i));
        
        %computes bi
        if q_table(i) == 1
            bi = min(a)+delta;  %arithmetic precision fix
        else
            bi = min(a)+delta:delta:max(a)-delta;
        end
        
        %compute yi
        yi = linspace(bi(1)-delta/2,bi(end)+delta/2,length(bi)+1);
        
        %        uncomment to optimize yi
        yi = zeros(1,numel(bi)+1);      %+-inf not included in bi
        for j = 1:length(yi)
            if j == 1
                yi(1) = sum(a(a<bi(j)))/sum(a<bi(j));
            elseif j == 2^q_table(i)
                yi(end) =  sum(a(a>bi(j-1)))/sum(bi(j-1)<a);
            else
                yi(j) = (-sum(a(a<bi(j-1)))+sum(a(a<bi(j))))/(-sum(a<bi(j-1))+sum(a<bi(j)));
            end
            if isnan(yi(j))
                yi(j) = bi(j)-delta;
            end
        end
        %        end of optimization%
        
        %quantization
        [~, c] = quantiz(a,bi,yi);
        quant((v_step*(mod(i-1,2^lev))+1:v_step*(mod(i-1,2^lev))+v_step),...
            (h_step*ceil(i/2^lev-1)+1:h_step*ceil(i/2^lev-1)+h_step)) = reshape(c,[v_step, h_step]);
        
        %entropy coding
        [p, sym] = hist(c,unique(c));   %find symbols and their probabilities (need to be normalized)
        dict = huffmandict(sym,p/numel(c));
        dictionary{i} = dict;
        encoded = cat(2,encoded,huffmanenco(c,dict));
        len(i) = length(encoded);
    end
end
end