function vec = fgradf(a0,v0,h0)

a = valder(a0); %angle in degrees
v = valder(v0); %velocity in ft/sec
h = valder(h0); %height in ft
rad = a*pi/180;
tana = tan(rad);
vhor = (v*cos(rad))^2;
f = (vhor/32)*(tana + sqrt(tana^2+64*h/vhor)); %In forward pass immediate partial derivatives are calculated and stored
mat = single(f); % we get immediate partial derivatives as well as a mapping for constructing extended Jacobian Matrix  

[m,n] = size(mat);
%Formation of extended Jacobian Matrix
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
%this is our extended Jacobian Matrix
new = new - eye(m);

derivatives = zeros(1,m);

%Reverse pass
derivatives(1,m) = 1;
  for j=m-1:-1:1
    for i=m:-1:j+1
	derivatives(1,j) = derivatives(1,j) + derivatives(1,i)*new(i,j);
    end
  end
f = double(f);
f(:,1)
%Finally we get the derivatives
derivatives(:,[1:3])
