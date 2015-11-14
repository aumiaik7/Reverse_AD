classdef valder2
   % valder2 class implementing Automatic Differentiation by operator overloading.
   % Computes first order derivative or multivariable gradient vectors by
   % starting with a known simple valder2 such as x=valder2(3,1) and 
   % propagating it through elementary functions and operators. 
   % by Richard D. Neidinger 10/23/08


   properties
      val  %function value
      der  %derivative value or gradient vector

     
   end 
   methods

      function obj = valder2(a,b)


          sprintf('**  valder2   **');
          a
          b
	  	
         %valder2 class constructor; only the bottom case is needed.
         if nargin == 0 %never intended for use.
            obj.val = [];
            #obj.der = [];
	    #obj.store = [];
         elseif nargin == 1 %c=valder2(a) for constant w/ derivative 0.
            obj.val = a;
            #obj.der = 0;
	  # obj.store = 0;
         else
            obj.val = a; %given function value

	    obj.der = [b];  %given derivative value or gradient vector
	    #obj.store = c;
         end
      end
      function vec = double(obj)
         %valder2/DOUBLE Convert valder2 object to vector of doubles.
         sprintf('**  double   **')
         vec = [ obj.der ];
      end
      function h = plus(u,v)
          sprintf('**  plus   **');
         %valder2/PLUS overloads addition + with at least one valder2 object argument
         #{
	 if ~isa(u,'valder2') %u is a scalar
            #h = valder2(u+v.val, v.der);
            h = valder2(u+v.val, 1);
         elseif ~isa(v,'valder2') %v is a scalar
            #h = valder2(u.val+v, u.der);
            h = valder2(u.val+v, 1);
         else
	#}
            #h = valder2(u.val+v.val, u.der+v.der);
            h = valder2(u.val+v.val, 2);
        # end
      end
      function h = uminus(u)
          sprintf('**  uminus   **');
         %valder2/UMINUS overloads negation - with a valder2 object argument
         #h = valder2(-u.val, -u.der);
         h = valder2(-u.val, -1);
      end
      function h = minus(u,v)
          sprintf('**  minus   **');
         %valder2/MINUS overloads subtraction - with at least one valder2 object argument
         if ~isa(u,'valder2') %u is a scalar
            #h = valder2(u-v.val, -v.der);
            h = valder2(u-v.val, -1);
         elseif ~isa(v,'valder2') %v is a scalar
            #h = valder2(u.val-v, u.der);
            h = valder2(u.val-v, 1);
         else
            h = valder2(u.val-v.val, 0);
         end
      end
      function h = mtimes(u,v)
          sprintf('**  mtimes   **')

         %valder2/MTIMES overloads multiplication * with at least one valder2 object argument
         if ~isa(u,'valder2') %u is a scalar
            #h = valder2(u*v.val, u*v.der);
            h = valder2(u*v.val, u, u);
             sprintf('**  1   **')
         elseif ~isa(v,'valder2') %v is a scalar
            #h = valder2(v*u.val, v*u.der);
            h = valder2(v*u.val, v ,v );
            sprintf('**  2   **')
         else
	   
            #h = valder2(u.val*v.val, u.der*v.val + u.val*v.der);
            h = valder2(u.val*v.val, [u.der v.der v.val u.val]);
            sprintf('**  3   **')
         end
      end
      function h = mrdivide(u,v)
         %valder2/MRDIVIDE overloads division / with at least one valder2 object argument
         sprintf('**  mdivide   **');
         if ~isa(u,'valder2') %u is a scalar
            #h = valder2(u/v.val, (-u*v.der)/(v.val)^2);
            h = valder2(u/v.val, (-u)/(v.val)^2);
         elseif ~isa(v,'valder2') %v is a scalar
            #h = valder2(u.val/v, u.der/v);
            h = valder2(u.val/v, 1/v);
         else
            #h = valder2(u.val/v.val, (u.der*v.val-u.val*v.der)/(v.val)^2);
            h = valder2(u.val/v.val, (v.val-u.val)/(v.val)^2);
         end
      end
      function h = mpower(u,v)
          sprintf('**  mpower   **');
         %valder2/MPOWER overloads power ^ with at least one valder2 object argument
         if ~isa(u,'valder2') %u is a scalar
            #h = valder2(u^v.val, u^v.val*log(u)*v.der);
            h = valder2(u^v.val, u^v.val*log(u));
         elseif ~isa(v,'valder2') %v is a scalar
            #h = valder2(u.val^v, v*u.val^(v-1)*u.der);
            h = valder2(u.val^v, v*u.val^(v-1));
         else
            h = exp(v*log(u)); %call overloaded log, * and exp
         end
      end
      function h = exp(u)
          sprintf('**  exp   **');
         %valder2/EXP overloads exp of a valder2 object argument
         #h = valder2(exp(u.val), exp(u.val)*u.der);
         h = valder2(exp(u.val), exp(u.val));
      end
      function h = log(u)
          sprintf('**  log   **');
         %valder2/LOG overloads natural logarithm of a valder2 object argument
         #h = valder2(log(u.val), (1/u.val)*u.der);
         h = valder2(log(u.val), (1/u.val));
      end
      function h = sqrt(u)
          sprintf('**  sqrt   **');
         %valder2/SQRT overloads square root of a valder2 object argument
         #h = valder2(sqrt(u.val), u.der/(2*sqrt(u.val)));
         h = valder2(sqrt(u.val), 1/(2*sqrt(u.val)));
      end
      function h = sin(u)
          sprintf('**  sin   **');
         %valder2/SIN overloads sine with a valder2 object argument
         #h = valder2(sin(u.val), cos(u.val)*u.der);
         h = valder2(sin(u.val), [u.der cos(u.val)]);
      end
      function h = cos(u)
          sprintf('**  cos   **');
         %valder2/COS overloads cosine of a valder2 object argument
         #h = valder2(cos(u.val), -sin(u.val)*u.der);
         h = valder2(cos(u.val), [u.der -sin(u.val)]);
      end
      function h = tan(u)
          sprintf('**  tan   **');
         %valder2/TAN overloads tangent of a valder2 object argument
         #h = valder2(tan(u.val), (sec(u.val))^2*u.der);
         h = valder2(tan(u.val), (sec(u.val))^2);
      end
      function h = asin(u)
          sprintf('**  asin   **');
         %valder2/ASIN overloads arcsine of a valder2 object argument
         #h = valder2(asin(u.val), u.der/sqrt(1-u.val^2));
         h = valder2(asin(u.val), 1/sqrt(1-u.val^2));
      end
      function h = atan(u)
          sprintf('**  atan   **');
         %valder2/ATAN overloads arctangent of a valder2 object argument
         #h = valder2(atan(u.val), u.der/(1+u.val^2));
         h = valder2(atan(u.val), 1/(1+u.val^2));
      end
   end
end
