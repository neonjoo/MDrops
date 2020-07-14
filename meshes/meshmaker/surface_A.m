function R = surface_A(b,lambda)
% https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.2016.0401

    a = lambda * b;
    
    f1 = @(x,y) sqrt((a^2*b^4 - (b^4 + 4*a^2*x.^2).^(3/2) + 4*a^2*x.^2.*sqrt(b^4 + 4*a^2*x.^2)) ...
             ./ (a^2*b^4 - (b^4 + 4*a^2*x.^2).^(3/2) + 4*a^2*x.^4 + 4*a^4*x.^2 + b^4*x.^2 + b^4*y.^2 + 4*a^2*x.^2.*y.^2));
    
    
    ymax = @(x) sqrt(sqrt(b^4+4*a^2*x.^2)-a^2-x.^2);
    xmax = sqrt(a^2+b^2);
    
    R = 1 - 8*quad2d(f1,0,xmax, 0, ymax );
end