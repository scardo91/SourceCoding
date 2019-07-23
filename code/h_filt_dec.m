function [out] = h_filt_dec(im, lev, type)
%computes the analisys stage of wavelet decomposition

%analisys row-wise
[L, H] = h_filt(im);

%analisys column-wise
[LL, LH] = h_filt(rot90(L));
[HL, HH] = h_filt(rot90(H));

%adjust results
LL = rot90(LL,3);
LH = rot90(LH,3);
HL = rot90(HL,3);
HH = rot90(HH,3);

%recursively apply the procedure for further decomposition
if lev > 1
    if strcmp(type,'uniform')
        LL = h_filt_dec(LL,lev-1,'uniform');
        LH = h_filt_dec(LH,lev-1,'uniform');
        HL = h_filt_dec(HL,lev-1,'uniform');
        HH = h_filt_dec(HH,lev-1,'uniform');
    elseif strcmp(type,'pyramidal')
        LL = h_filt_dec(LL,lev-1,'pyramidal');
    else
        error('unsupported decompostion strategy');
    end
    
end
%%%%uncomment only to adjust visualization of decomposed images
% LL = 255/max(LL(:)).*LL;
% HL = 255/max(HL(:)).*HL;
% LH = 255/max(LH(:)).*LH;
% HH = 255/max(HH(:)).*HH;

out = [LL, LH; HL, HH];
end