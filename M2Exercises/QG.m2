p = 3
F = GF(p);
R = F[x, y, z];
I1 = (ideal(x,y,z))^3;
monoms = gens I1

randomPoly = () -> (
    L = toList apply(1..10, i -> random(F));
    coeffs = transpose(matrix{L});
    return (monoms * coeffs)_(0,0);
);

a = 0
s = 10

for i in (1..s) do (
f = randomPoly();
if (f^2 % ideal(x^3, y^3, z^3) == 0) then (a = a +1);
);

print (a / s);


loadPackage "TestIdeals"
loadPackage "FrobeniusThresholds"

p = 7;
e = 3;
F = GF(p);
R = F[x,y];
I = ideal(x^2+y^3);
J = I;

s = 3
for i in (1..p^e) do (
	if (frobeniusRoot(e, J) != frobeniusRoot(e, I*J)) then (
		print(i);
	);
	J = I*J;
);


frobeniusNu(e, I, ideal(x,y))



load "NuInvariants.m2";
loadPackage "TestIdeals";


p = 7;
e = 2;
F = GF(p);
R = F[x,y];
I = ideal(x^2+y^3, y);

L = nuInvariants(p,e,I);
