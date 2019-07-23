function [reco] = h_filt_synt(deco, lev, type)
%reconstructs the decomposed image

[h, w] = size(deco);

if strcmp(type, 'uniform')
    while lev > 0
        %extract 4x4 subband's subsets
        hor = ones(1,2^(lev-1))*w/2^(lev-1);
        vert = ones(1,2^(lev-1))*h/2^(lev-1);
        sub = mat2cell(deco, vert, hor);
        
        %merge 4x4 block
        for i = 1:4^(lev-1)
            h_step = w/2^(lev-1); %hor size of the image ricomposed by 4 subbands
            v_step = h/2^(lev-1);
            ns = 2^(lev-1); %number of squares in a column/row
            
            %sinthesys column-wise
            L = h_filt_s(sub{i}(1:h/2^lev,1:h_step/2),sub{i}(1:h/2^lev,h_step/2+1:end));
            H = h_filt_s(sub{i}(h/2^lev+1:end,1:h_step/2),sub{i}(h/2^lev+1:end,h_step/2+1:end));
            
            %sinthesys row-wise
            deco((v_step*(mod(i-1,ns))+1:v_step*(mod(i-1,ns))+v_step),...
                (h_step*ceil(i/ns-1)+1:h_step*ceil(i/ns-1)+h_step)) = rot90(h_filt_s(rot90(L,3),rot90(H,3)));
        end
        %repeat untill only one 4x4 block
        reco = deco;
        lev = lev-1;
    end
    
elseif strcmp(type,'pyramidal')
    l = 1;
    reco = deco;
    
    %perform uniform reconstruction from the smallest (upper left) 4x4
    %block untill full image
    while ~isequal(l, lev+1)
        reco(1:h/(2^(lev-l)),1:w/(2^(lev-l))) = h_filt_synt(reco(1:h/(2^(lev-l)),1:w/(2^(lev-l))), 1, 'uniform');
        l = l+1;
    end
else
    error('unsupported decompostion strategy');
end
end