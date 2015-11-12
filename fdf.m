function vec = fdf(a)
%FDF takes a scalar and returns the double vector [ f(a), f'(a) ]
%   where f is defined in normal syntax below.
x = valder(a,1);
y = exp(-sqrt(x))*sin(x*log(1+x^2));
vec = double(y);