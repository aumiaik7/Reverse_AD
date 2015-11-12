function vec = fgradf(a0,v0,h0)
%FGRADF takes 3 scalars and returns the double vector 
%   [ f(a,b,c), df/dx(a,b,c), df/dy(a,b,c), df/dz(a,b,c) ]
%   where f is defined in normal syntax below.
%   This example is range of a tennis serve.
#{
a = valder(a0,[1 0 0]); %angle in degrees
v = valder(v0,[0 1 0]); %velocity in ft/sec
h = valder(h0,[0 0 1]); %height in ft
rad = a*pi/180;
tana = tan(rad);
vhor = (v*cos(rad))^2;
f = (vhor/32)*(tana + sqrt(tana^2+64*h/vhor)); %horizontal range
vec = double(f);
#}

x = valder(a0,[1]); %angle in degrees
#x = a0;
y = valder(v0,[1]); %velocity in ft/sec
#y = v0;
#z = h0;
z = valder(h0,[1]); %height in ft
#x = x*pi/180;
#y = y*pi/180;
f = x*z*sin(x*y);

vec = double(f);
