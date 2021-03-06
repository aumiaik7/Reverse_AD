classdef valder
   % VALDER class implementing Automatic Differentiation by operator overloading.
   % Computes first order derivative or multivariable gradient vectors by
   % starting with a known simple valder such as x=valder(3,1) and 
   % propagating it through elementary functions and operators. 
   % by Richard D. Neidinger 10/23/08
	
   properties
      val  %function value
      der  %derivative value or gradient vector
      pos  %variable position in forward accumulation
      
   end 
   methods

      function obj = valder(a,b)
          sprintf('**  valder   **')
	 
   	  global store = [];
	
	%global ar
         %VALDER class constructor; only the bottom case is needed.
         if nargin == 0 %never intended for use.
            obj.val = [];
     
   
         elseif nargin == 1 %c=valder(a) for constant w/ derivative 0.
            obj.val = a;
	    [row,col] = size(store);
	    store = [store; 0 0 row 0]; %stores the value and position for forward accumulation
	    obj.pos = row + 1; % position
	    
         else
	    obj.val = a; %given function value
            [row,col] = size(store);
	    store = [store;b]; %stores the value and position for forward accumulation
	    obj.pos = row + 1;% position
	 end
      end

      function mat = single(obj)
	global store;
	mat = store;
      end	
      function vec = double(obj)
         %VALDER/DOUBLE Convert valder object to vector of doubles.
         sprintf('**  double   **') ; 
         vec = [ obj.val obj.der ]
      end
      function h = plus(u,v)
          sprintf('**  plus   **')
         %VALDER/PLUS overloads addition + with at least one valder object argument
         if ~isa(u,'valder') %u is a scalar
            %h = valder(u+v.val, v.der)
	    h = valder(u+v.val, [0 1 0 v.pos]);
         elseif ~isa(v,'valder') %v is a scalar
            % h = valder(u.val+v, u.der)
	    h = valder(u.val+v, [1 0 u.pos 0]);
         else
            %h = valder(u.val+v.val, u.der+v.der)
	    h = valder(u.val+v.val, [1 1 u.pos v.pos]);
         end
      end
      function h = uminus(u)
          sprintf('**  uminus   **');
         %VALDER/UMINUS overloads negation - with a valder object argument
        % h = valder(-u.val, -u.der)
	   h = valder(-u.val, [-1 0 u.pos 0]);
      end
      function h = minus(u,v)
          sprintf('**  minus   **');
         %VALDER/MINUS overloads subtraction - with at least one valder object argument
         if ~isa(u,'valder') %u is a scalar
            %h = valder(u-v.val, -v.der)
	    h = valder(u-v.val, [0 -1 0 v.pos]);
         elseif ~isa(v,'valder') %v is a scalar
           % h = valder(u.val-v, u.der)
	    h = valder(u.val-v, [1 0 u.pos 0]);
         else
            %h = valder(u.val-v.val, u.der-v.der)
		h = valder(u.val-v.val, [1 -1 u.pos -v.pos]);
         end
      end
     function h = mtimes(u,v)
          sprintf('**  mtimes   **');

         %VALDER/MTIMES overloads multiplication * with at least one valder object argument
         if ~isa(u,'valder') %u is a scalar
            #h = valder(u*v.val, u*v.der);
            h = valder(u*v.val, [0 u 0 v.pos]);
             
         elseif ~isa(v,'valder') %v is a scalar
            #h = valder(v*u.val, v*u.der);
            h = valder(v*u.val, [v 0 u.pos 0]);
         else
	   
            #h = valder(u.val*v.val, u.der*v.val + u.val*v.der);
            h = valder(u.val*v.val, [v.val u.val u.pos v.pos]);
            sprintf('**  3   **')
         end
      end
      function h = mrdivide(u,v)
         %VALDER/MRDIVIDE overloads division / with at least one valder object argument
         sprintf('**  mdivide   **');
         if ~isa(u,'valder') %u is a scalar
            %h = valder(u/v.val, (-u*v.der)/(v.val)^2)
	    h = valder(u/v.val, [ 0 (-u)/(v.val)^2 0 v.pos]);
         elseif ~isa(v,'valder') %v is a scalar
            %h = valder(u.val/v, u.der/v)
	    h = valder(u.val/v, [1/v 0 u.pos 0]);
         else
            %h = valder(u.val/v.val, (u.der*v.val-u.val*v.der)/(v.val)^2)
            h = valder(u.val/v.val, [1/v.val (-u.val)/(v.val)^2 u.pos v.pos]);
         end
      end
      function h = mpower(u,v)
          sprintf('**  mpower   **');
         %VALDER/MPOWER overloads power ^ with at least one valder object argument
         if ~isa(u,'valder') %u is a scalar
            %h = valder(u^v.val, u^v.val*log(u)*v.der)
	    h = valder(u^v.val, [0 u^v.val*log(u) 0 v.pos]);
         elseif ~isa(v,'valder') %v is a scalar
            %h = valder(u.val^v, v*u.val^(v-1)*u.der)
	    h = valder(u.val^v, [v*u.val^(v-1) 0 u.pos 0]);
         else
            h = exp(v*log(u)); %call overloaded log, * and exp
         end
      end
      function h = exp(u)
          sprintf('**  exp   **');
         %VALDER/EXP overloads exp of a valder object argument
         % h = valder(exp(u.val), exp(u.val)*u.der)
	 h = valder(exp(u.val), [exp(u.val) 0 u.pos 0]);
      end
      function h = log(u)
          sprintf('**  log   **');
         %VALDER/LOG overloads natural logarithm of a valder object argument
         %h = valder(log(u.val), (1/u.val)*u.der)
 	 h = valder(log(u.val), [1/u.val 0 u.pos 0]);
      end
      function h = sqrt(u)
          sprintf('**  sqrt   **');
         %VALDER/SQRT overloads square root of a valder object argument
         %h = valder(sqrt(u.val), u.der/(2*sqrt(u.val)))
	 h = valder(sqrt(u.val), [1/(2*sqrt(u.val)) 0 u.pos 0]);
      end
      function h = sin(u)
          sprintf('**  sin   **');
         %VALDER/SIN overloads sine with a valder object argument
         %h = valder(sin(u.val), cos(u.val)*u.der)
	h = valder(sin(u.val), [cos(u.val) 0 u.pos 0]);
      end
      function h = cos(u)
          sprintf('**  cos   **');
         %VALDER/COS overloads cosine of a valder object argument
         %h = valder(cos(u.val), -sin(u.val)*u.der)
         h = valder(cos(u.val), [-sin(u.val) 0 u.pos 0]);
      end
      function h = tan(u)
          sprintf('**  tan   **');
         %VALDER/TAN overloads tangent of a valder object argument
         %h = valder(tan(u.val), (sec(u.val))^2*u.der)
	h = valder(tan(u.val), [(sec(u.val))^2 0 u.pos 0]);
      end
      function h = asin(u)
          sprintf('**  asin   **');
         %VALDER/ASIN overloads arcsine of a valder object argument
         %h = valder(asin(u.val), u.der/sqrt(1-u.val^2))
         h = valder(asin(u.val), [1/sqrt(1-u.val^2) 0 u.pos 0]);
      end
      function h = atan(u)
          sprintf('**  atan   **');
         %VALDER/ATAN overloads arctangent of a valder object argument
         %h = valder(atan(u.val), u.der/(1+u.val^2))
         h = valder(atan(u.val), [1/(1+u.val^2) 0 u.pos 0]);
      end
   end
end
