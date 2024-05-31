function para = parameters(mat)

% Compute the 8 parameters based on the association matrix (for a specific
% patient, indicator and bWidth)

% 1.Clustering coefficient 2.Characteristic path length 3.Strength 4.Global efficiency   

mat = squeeze(mat);
mat = mat+mat';
for i = 1:size(mat,1)
    mat(i,i) = 0;
end
normMat = weight_conversion(abs(mat), 'normalize');

para{1} = clustering_coef_wu(mat);
para{2} = charpath(distance_wei(1./mat),0,0);
para{3} = strengths_und(mat);
para{4} = efficiency_wei(mat,0);