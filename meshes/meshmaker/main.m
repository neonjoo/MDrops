% a = 0.9;
% b = 0.9;
% c = 1/(a*b);
%fd=@(p)  p(:,1).^2/a^2+p(:,2).^2/b^2+p(:,3).^2/c^2-1;
lambda = 0.8;
bfun = @(bvar) surface_A(bvar,lambda);
b = fzero(bfun,1);
a = lambda * b;

fd=@(p)  (p(:,1).^2 + p(:,2).^2 + p(:,3).^2).^2 ...
        - 2*a^2*(p(:,1).^2 - p(:,2).^2 - p(:,3).^2) ...
        + a^4 - b^4;

triangsize = 0.15;

xlim = 2.6;
ylim = 2.6;
zlim = 2.6;
%p - node positions (3N); t - triangle indicies (3N)
[p,t]=distmeshsurface(fd,@huniform,triangsize,[-xlim,-ylim,-zlim; xlim,ylim,zlim]);
%%
S = 0;

for i = 1:size(t,1) 
   x1 = p(t(i,1),:);
   x2 = p(t(i,2),:);
   x3 = p(t(i,3),:);
   
   S = S + 0.5 * norm(cross(x3-x1,x2-x1));
end

%csvwrite('points_ellpise_fewN.csv',p)
%csvwrite('faces_ellipse_fewN.csv',t)