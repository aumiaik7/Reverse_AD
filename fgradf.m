function vec = fgradf(a0,v0,h0)
%FGRADF takes 3 scalars and returns the double vector 
%   [ f(a,b,c), df/dx(a,b,c), df/dy(a,b,c), df/dz(a,b,c) ]
%   where f is defined in normal syntax below.
%   This example is range of a tennis serve.

a = valder(a0); %angle in degrees
v = valder(v0); %velocity in ft/sec
h = valder(h0); %height in ft
rad = a*pi/180;
tana = tan(rad);
vhor = (v*cos(rad))^2;
f = (vhor/32)*(tana + sqrt(tana^2+64*h/vhor)); %horizontal range
#vec = double(f);
mat = single(f);
#}
#{
x = valder(a0); 
y = valder(v0); 
z = valder(h0); 
v1 = x*z;
v2 = x*y;
v3 = sin(v2);
v4 = v1*v3;
mat = single(v4)
#}

[m,n] = size(mat);
new = zeros(m,m);

  for i=1:m
	
    if mat(i,1) ~= 0	 
	j = int16(mat(i,3));
    	new(i,j) = mat(i,1);
    end
    if mat(i,2) ~= 0
	j = int16(mat(i,4));	 
    	new(i,j) = mat(i,2);
    end 
  end
new = new - eye(m);
derivatives = zeros(1,m);
derivatives(1,m) = 1
  for j=m-1:-1:1
    for i=m:-1:j+1
	
	derivatives(1,j) = derivatives(1,j) + derivatives(1,i)*new(i,j);
    end
  end

derivatives
