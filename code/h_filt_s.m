function [rec] = h_filt_s(to_mer1,to_mer2)
%computes column-wise recontruction

%compute recostruction filter coefficients
[~,~,ILp,IHp] = wfilters('haar');

%initialize quantities
[h, w] = size(to_mer1);
L = numel(ILp);
rec = zeros(2*h+3*L-1,w);
high = zeros(h+L,w);
low =  high;

%pad and upsample along columns
for i = 1:w
    high(:,i) = padarray(to_mer2(:,i),L/2,'circular');
    low(:,i) = padarray(to_mer1(:,i),L/2,'circular');
end
high = upsample(high,2);
low = upsample(low,2);

%filter
for i = 1:w
    rec(:,i) = conv(double(low(:,i)),ILp) + conv(double(high(:,i)),IHp);
end

%select the right values
rec = rec(length(IHp)+1:2*h+length(IHp),:);
end