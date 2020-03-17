function [ outs ] = null_spec_polynomial( G, r, z )
%Evaluate the generalized null spectrum function of given G. And it's not a
%polynomial anymore

Z_mat       = z*ones(size(G));          %matrix holding z values
EDM         = (squareform(pdist(r(:)))).*(-1*ones(size(G))+triu(2*ones(size(G))));    %correct sign EDM
outs        = abs(sum( G(:).*(Z_mat(:).^EDM(:)))); %null spectrum polynomial at z


end

