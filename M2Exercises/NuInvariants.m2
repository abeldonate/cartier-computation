loadPackage "TestIdeals";


nuInvariants = method();

nuInvariants (ZZ, ZZ, Ideal) := (p, e, I) -> (
	s := numgens I;
        nu := ();
	J := I;
        for i in (1..s*p^e) do (
		if (frobeniusRoot(e, J) != frobeniusRoot(e, I*J)) then ( 
			nu = append(nu, i);
		);
		J = I*J;
	);
	return nu;
);



p = 7;
e = 2;
F = GF(p);
R = F[x, y];
I = ideal(x^2+y^3, y);

L = nuInvariants(7,2,I);
print L;

