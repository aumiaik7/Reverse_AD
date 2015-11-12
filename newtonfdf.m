function root = newtonfdf(a)
%NEWTON seeks a zero of the function defined in fdf using the initial a
% root estimate and Newton's method (with no exception protections).
% fdf uses @valder to return a vector of function and derivative values.
delta = 1;
while abs(delta) > .000001
    fvec = fdf(a);
    delta = fvec(1)/fvec(2); %value/derivative
    a = a - delta;
end
root = a;