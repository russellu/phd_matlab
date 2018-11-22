% Alternative way of PCA computaion for large feature k and small
% obervation
function [COEFF1, SCORE1, latent1]=PCA_alter(X)
% X is the input matrix in the dimension n*k (observation by dimension)
X_center=X-repmat(mean(X), size(X,1),1);
[S V D]=svd(X_center*X_center');
COEFF1=X_center'*S./repmat(sqrt(diag(V))',size(X,2),1);
SCORE1=X_center*COEFF1;
latent1=diag(cov(SCORE1));