function stacked =  vecangle(complvect)
%VECANGLE Takes a complex vector of size n and stacks real and imaginary part to result in a vector if size 2n.
stacked = [real(complvect), imag(complvect)];
end

