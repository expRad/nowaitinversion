function [f3,d3, lamb, mu] = line_step(f1, f2, d1, d2, damping, compl)
% performs the linear step, while treating the input as a single large
% vector, which defines the linear space

% 
% if ~exist(compl)
%     compl = true;
% end
% 
if (norm(f2(:)) < 1e-10) || (norm(d2(:)) < 1e-10)
    error("f2==0 or d2==0, therefore no linear space could be defined")
end


asize = size(squeeze(f1));
a = flatten(f1);
x = flatten(f2);

bsize = size(squeeze(d1));
b = flatten(d1);
y = flatten(d2);

if compl
    Rya = dot(real(y), real(a));
    Ryb = dot(real(y), real(b));
    Ryy = dot(real(y), real(y));
    Rxa = dot(real(x), real(a));
    Rxb = dot(real(x), real(b));
    Rxx = dot(real(x), real(x));
    Rxy = dot(real(x), real(y));
    
    Iya = dot(imag(y), imag(a));
    Iyb = dot(imag(y), imag(b));
    Iyy = dot(imag(y), imag(y));
    Ixa = dot(imag(x), imag(a));
    Ixb = dot(imag(x), imag(b));
    Ixx = dot(imag(x), imag(x));
    Ixy = dot(imag(x), imag(y));
    
    ABY = Rya + Iya - Ryb - Iyb;
    ABX = Rxa + Ixa - Rxb - Ixb;
    XX = Rxx + Ixx;
    YY = Ryy + Iyy;
    XY = Rxy + Ixy;
    

    lamb = ( - YY*ABX + ABY.*XY) / (YY * XX - XY.^2); 
    mu = (ABY + lamb*XY) / YY;
else
    % for real a,b,x,y!!
    
    a = stackcompl(a);
    b = stackcompl(b);
    x = stackcompl(x);
    y = stackcompl(y);
        
    ABY = dot(a,y) - dot(b,y);
    BAX = dot(b,x) - dot(a,x);
    SX = dot(x,x);
    SY = dot(y,y);
    XY = dot(x,y);
%     lamb = (ABY.*XY / SY + BAX) / (SX - XY.^2 / SY);
    lamb = (ABY.*XY + SY*BAX) / (SY*SX - XY.^2 );
    mu = lamb*(XY/SY) + ABY / SY;
    
    a = unstackcompl(a);
    b = unstackcompl(b);
    x = unstackcompl(x);
    y = unstackcompl(y);
    
end

f3 =  a + damping*lamb*x;
d3 = b + damping*mu*y;

f3 = reshape(f3, asize);
d3 = reshape(d3, bsize);

% if isreal(f1) && isreal(f2) && isreal(d1) && isreal(d2)
% else

end