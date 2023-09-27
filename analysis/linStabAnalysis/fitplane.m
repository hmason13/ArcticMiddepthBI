function [p,N] = fitplane(x,y,f)

% [p,N] = fitplane(x,y,f) Find best fit plane p (with normal N) to
% surface f = f(x,y), with size(f) = [length(x),length(y)].  Output
% p has same dimensions as f.  Normal N is a 3-vector such that 
% N.X = 1. (so slope_x -N(1)/N(3) and slope_y = -N(2)/N(3))

if (size(f,1))~=length(x), error('size(x) neq size(f,1)'), end
if (size(f,2))~=length(y), error('size(y) neq size(f,2)'), end

[x_,y_] = ndgrid(x,y);
nt = numel(x_);  
xl = reshape(x_,[nt 1]);
yl = reshape(y_,[nt 1]);
zl = reshape(f ,[nt 1]);
R = [xl(:) yl(:) zl(:)];
r0 = 2*min(R)-max(R)-1;
N = (R-r0(ones(nt,1),:))\ones(nt,1);
N = N/(1+r0*N);
%[N,d] = fitplane(xl,yl,zl);
p = (1/N(3)) - (N(1)/N(3))*x_ - (N(2)/N(3))*y_;
