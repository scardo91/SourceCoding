function map = bit_alloc_unif(img, R, lev)
% works only for uniform decomposition and computes the optimal bit
% allocation using iterative procedure

[h, w] = size(img);
M = 4^lev;  %number of subbands
Rb = M*R;     %total available bits
map = zeros(1,M);

%extract subbands
hor = ones(1,2^lev)*w/2^lev;
vert = ones(1,2^lev)*h/2^lev;
sub = mat2cell(img, vert, hor);

%initialization
sigma_rk = zeros(1,M);
for i = 1:M
    sigma_rk(i) = 2^(-2*map(1,i))*var(sub{i}(:));
end

%iterative allocation
while Rb ~= 0
    [~ , ind] = sort(sigma_rk,'descend');
    map(1,ind(1)) = map(1,ind(1)) +1;
    sigma_rk(ind(1)) = sigma_rk(ind(1))/4;
    Rb = Rb - 1;
end

%alternative strategy
% sigma_rk = zeros(1,M);
% geom_mean = 1;
% for i = 1:M
%     geom_mean = geom_mean*var(sub{i}(:));
% end
% geom_mean=geom_mean^(1/M);
% for j = 1:M
%     map(j) = max(0,round(R + (1/2)*log2(var(sub{j}(:))/geom_mean)));
% end
% if sum(map) > M
%     a = find(map == 1,sum(map)-M,'last');
%     map(a) = 0;
% end

end