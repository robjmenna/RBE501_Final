function H = generate_H(a, alpha, d, theta)

H = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
iterations = size(a,2);
for i = 1:iterations
H = H*DH(alpha(i), a(i), d(i), theta(i));
end
end