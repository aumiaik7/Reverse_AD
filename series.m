classdef series
   % SERIES class implementing AD to compute series coefficients.
   % Start with x = series(a,1,d) for variable x at value a through degree d.
   % Propagate forward through elementary operations and functions overloaded 
   % in the methods below.
   % by Richard Neidinger 11/17/08
   properties
      val  %function value (constant term)
      coef %vector of Taylor coefficients, linear to highest term
   end 
   methods
      function obj = series(a,der,order)
         %SERIES series class constructor
         %   f = series(a,der,order) creates a series object for a Talyor polynomial with 
         %   constant term a, following coefficient(s) in der, and zeros thru degree order.
         %   h = series(a,der) will assume that der contains all desired coefficients.
         %   x = series(a,1,order) creates a series for variable x at value a.
         if nargin == 0 % never intended for use.
            obj.val = [];
            obj.coef = [];
         elseif isa(a,'series')  % pass through if already a series object.
             obj = a;
         elseif nargin == 1  % never intended for use.
             obj.val = a;
             obj.coef = [];
         else
             obj.val = a;
             obj.coef = der;
             if nargin >= 3 
                 obj.coef = [obj.coef, zeros(1,order-length(der))];
             end
         end
      end
      function vector = double(obj)
         %SERIES/DOUBLE Convert series object to vector of doubles
         %   containing all the series coefficient values.
         vector = [ obj.val, obj.coef ];
      end
      function h = plus(u,v)
         %SERIES/PLUS overloads addition with at least one series object argument
         % assumes .coef lists are of same length
         if ~isa(u,'series')
             h = series(u + v.val, v.coef);
         elseif ~isa(v,'series')
             h = series(v + u.val, u.coef);
         else
             h = series(u.val + v.val, u.coef + v.coef);
         end
      end
      function h = uminus(u)
         %SERIES/UMINUS overloads negation with a series object argument
         h = series(-u.val, -u.coef);
      end
      function h = minus(u,v)
         %SERIES/MINUS overloads subtraction with at least one series object argument
         % assumes .coef lists are of same length
         if ~isa(u,'series')
             h = series(u - v.val, -v.coef);
         elseif ~isa(v,'series')
             h = series(u.val - v, u.coef);
         else
             h = series(u.val - v.val, u.coef - v.coef);
         end
      end
      function h = mtimes(u,v)
         %SERIES/MTIMES overloads multiplication with at least one series object argument
         % assumes .coef lists are of same length
         if ~isa(u,'series') % u is a scalar
             h = series(u*v.val, u*v.coef);
         elseif ~isa(v,'series')  % v is a scalar
             h = series(v*u.val, v*u.coef);
         else
             d = length(u.coef);
             uvec = [u.val, u.coef];
             vvec = [v.val, v.coef];
             h = series(u.val*v.val,0,d);
             for k = 1:d
                 h.coef(k) = uvec(1:k+1) * vvec(k+1:-1:1)' ;
             end
         end
      end
      function h = mrdivide(u,v)
         %SERIES/MRDIVIDE overloads division with at least one series object argument
         % assumes .coef lists are of same length
         if ~isa(v,'series')  % assume v is a scalar
             h = series(u.val/v, u.coef/v);
         else
             d = length(v.coef);
             if ~isa(u,'series')  % assume u is a scalar
                 u = series(u,0,d);
             end
             hvec = zeros(1,d+1);  % hvec(k+1) to hold h_k series coef
             hvec(1) = u.val/v.val;
             for k = 1:d
                 hvec(k+1) = (u.coef(k) - hvec(1:k) * v.coef(k:-1:1)') /v.val;
             end
             h = series(hvec(1), hvec(2:d+1));
         end
      end
      function h = sqrt(u)
         %SERIES/SQRT overloads square root of a series object argument
         d = length(u.coef);
         h = series(sqrt(u.val),0,d);
         h02 = 2*h.val;
         for k = 1:d
             h.coef(k) = (u.coef(k) - h.coef(1:k-1) * h.coef(k-1:-1:1)') / h02;
         end
      end
      function h = exp(u)
         %SERIES/EXP overloads exp with a series object argument
         d = length(u.coef);
         up = (1:d).*u.coef;  %vector of coefficients for u'(x)
         hvec = [exp(u.val), zeros(1,d)]; %hvec(k+1) holds h_k coef
         for k = 1:d
             hvec(k+1) = (up(1:k) * hvec(k:-1:1)') / k;
         end
         h = series(hvec(1), hvec(2:d+1));
      end
      function h = log(u)
         %SERIES/LOG overloads natural logarithm of a series object argument
         d = length(u.coef);
         hp = [u.coef(1)/u.val, zeros(1,d-1)];  % kth entry to hold k times series coef h_k.
         for k = 2:d
             hp(k) = (k*u.coef(k) - hp(1:k-1) * u.coef(k-1:-1:1)') / u.val;
         end
         h = series( log(u.val), hp./(1:d) );
      end
      function h = mpower(u,r)
         %SERIES/MPOWER overloads power with at least one series object argument.
         % Assumes .coef lists are of same length, if both are series objects.
         if isa(r,'series')
             h = exp(r*log(u));  % if u is scalar, only exp is series "convolution"
                                 % if u and r are series, exp, log, and * are "convolutions"
         else % assume r is a scalar
             if r == 2
                 h = u*u;
             elseif r == 1/2
                 h = sqrt(u);
             else % about same computation as exp&log but this allows negative base
                 d = length(u.coef);
                 up = (1:d).*u.coef;
                 hvec = [u.val^r, zeros(1,d)];  % hvec(k+1) to hold h_k series coef
                 hp = zeros(1,d);  % kth entry to hold k times series coef h_k.
                 for k = 1:d
                     hp(k) = (r*up(1:k)*hvec(k:-1:1)' - hp(1:k-1)*u.coef(k-1:-1:1)') /u.val;
                     hvec(k+1) = hp(k)/k;
                 end
                 h = series(hvec(1), hvec(2:d+1));
             end
         end
      end
      function [s, c] = sincos(u)
         %SERIES/SINCOS overloads exp with a series object argument
         d = length(u.coef);
         up = (1:d).*u.coef;
         svec = [sin(u.val), zeros(1,d)];
         cvec = [cos(u.val), zeros(1,d)];
         for k = 1:d
             svec(k+1) = ( up(1:k) * cvec(k:-1:1)' ) / k;
             cvec(k+1) = -( up(1:k) * svec(k:-1:1)' ) / k;
         end
         s = series(svec(1), svec(2:d+1));
         c = series(cvec(1), cvec(2:d+1));
      end
      function h = sin(u)
         %SERIES/SIN overloads sine with a series object argument
         [h, g] = sincos(u);
      end
      function g = cos(u)
         %SERIES/COS overloads cosine of a series object argument
         [h, g] = sincos(u);
      end
      function h = tan(u)
         %SERIES/TAN overloads tangent of a series object argument
         d = length(u.coef);
         up = (1:d).*u.coef;
         hvec = [tan(u.val), zeros(1,d)];  % hvec(k+1) to hold h_k series coef
         hp = zeros(1,d);  % kth entry to hold k times series coef h_k.
         vvec = [(cos(u.val))^(-2),zeros(1,d)];  % auxiliary v= sec(u)^2= 1+h^2
         for k = 1:d
             hp(k) = up(1:k)* vvec(k:-1:1)';
             hvec(k+1) = hp(k)/k;
             vvec(k+1) = 2*(hp(1:k)*hvec(k:-1:1)')/k;
         end
         h = series(hvec(1), hvec(2:d+1));
      end
      function h = asin(u)
         %SERIES/ASIN overloads arcsine of a series object argument
         d = length(u.coef);
         uvec = [u.val, u.coef];
         hp = zeros(1,d);  % kth entry to hold k times series coef h_k.
         v0 = sqrt(1-u.val^2);  % auxiliary v= cos(h)= sqrt(1-u^2)
         vcoef = zeros(1,d);
         for k = 1:d
             hp(k) = (k*u.coef(k) - hp(1:k-1)*vcoef(k-1:-1:1)') /v0;
             vcoef(k) = -(hp(1:k)*uvec(k:-1:1)')/k;
         end
         h = series(asin(u.val), hp./(1:d));
      end
      function h = atan(u)
         %SERIES/ATAN overloads arctangent of a series object argument
         d = length(u.coef);
         uvec = [u.val, u.coef];
         up = (1:d).*u.coef;
         hp = zeros(1,d);  % kth entry to hold k times series coef h_k.
         v0 = 1+u.val^2;  % auxiliary v= 1+u^2 = sec(h)^2
         vcoef = zeros(1,d);
         for k = 1:d
             hp(k) = (up(k) - hp(1:k-1)*vcoef(k-1:-1:1)') /v0;
             vcoef(k) = 2*(up(1:k)*uvec(k:-1:1)')/k;
         end
         h = series(atan(u.val), hp./(1:d));
      end
   end
end