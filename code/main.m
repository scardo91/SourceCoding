clear all;
close all;

levels = 3;
strategy = 'uniform';
R = 1;
encoding = true;
original = imread('barbara.png');
[h, w] = size(original);

if strcmp(strategy,'pyramidal')
    encoding = false;
end

%decompose image
decomposed = h_filt_dec(original, levels, strategy);

%compute coding gain
Gain = comp_gain(decomposed, levels, strategy);
H = entropy(original); %measure of bits/pixel of the original image once
                        %entropy coded

% performs quantization on the original image (only for R integer)
if mod(R,2) == 1 || mod(R,2) == 0
    
    a = double(reshape(original,[1,numel(original)]));
    delta = (max(a)-min(a))/2^(R);
    bi = min(a)+delta:delta:max(a)-delta;
    yi = bi(1)-delta/2:delta:bi(end)+delta/2;
    
    %uncomment to optimize recostruction levels%
    yi = zeros(1,numel(bi)+1);      %+-inf not included in bi
    for j = 1:length(yi)
        if j == 1
            yi(1) = sum(a(a<bi(j)))/sum(a<bi(j));
        elseif j == 2^R
            yi(end) =  sum(a(a>bi(j-1)))/sum(bi(j-1)<a);
        else
            yi(j) = (-sum(a(a<bi(j-1)))+sum(a(a<bi(j))))/(-sum(a<bi(j-1))+sum(a<bi(j)));
        end
        if isnan(yi(j))
            yi(j) = bi(j)-delta;
        end
    end
    %end of yi optimization%
    
    [~, c] = quantiz(a,bi,yi);
    c = round(c);
    orig_quant = uint8(reshape(c,[h, w]));
    
    diff1 = (orig_quant-original);
    MSE = comp_mse(diff1);
end

%the following works only for uniform decomposition

%compute quantization table
q_table = bit_alloc_unif(decomposed,R,levels);

%quantize and entropy coding
[quant_decomposed, coded, dictionaries, lengths] = quant_subband(decomposed,q_table);

%ricompose the original image after quantization (comment to test entropy decoding)
ricomposed = h_filt_synt(quant_decomposed,levels,strategy);

%%uncomment to test reconstruction without quantizing
%ricomposed = h_filt_synt(decomposed,levels,strategy);

if encoding
    %entropy decoding
    M = 4^levels;
    v_step = h/2^(levels);
    h_step = w/2^(levels);
    decoded_quant = zeros(h,w);
    lengths = cat(2,0,lengths);
    for i = 1:M
        if ~iscell(dictionaries{i})
            decoded_quant((v_step*(mod(i-1,2^levels))+1:v_step*(mod(i-1,2^levels))+v_step),...
                (h_step*ceil(i/2^levels-1)+1:h_step*ceil(i/2^levels-1)+h_step)) = zeros(v_step,h_step);
        else
            decoded_quant((v_step*(mod(i-1,2^levels))+1:v_step*(mod(i-1,2^levels))+v_step),...
                (h_step*ceil(i/2^levels-1)+1:h_step*ceil(i/2^levels-1)+h_step)) = ...
                reshape(huffmandeco(coded(lengths(i)+1:lengths(i+1)),dictionaries{i}),[v_step,h_step]);
        end
    end
    ricomposed = h_filt_synt(decoded_quant,levels,strategy);
end

%figures
% iptsetpref('ImshowBorder','tight');
ricomposed = uint8(ricomposed);
if mod(R,2) == 1 || mod(R,2) == 0
    subplot(2,2,1)
    imshow(original);
    title('Original Image.')
    subplot(2,2,2)
    imshow(orig_quant);
    title('Direct Quantization.')
    subplot(2,2,3)
    imshow(ricomposed);
    title('SBC Quantization.')
else
    subplot(1,2,1)
    imshow(original);
    title('Original Image.')
    subplot(1,2,2)
    imshow(ricomposed);
    title('SBC Quantization.')
end


%compute distorsion and real gain (if possible)
diff2 = ricomposed-original;
MSE2 = comp_mse(diff2);
BPP = numel(coded)/(w*h);