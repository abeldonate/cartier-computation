loadPackage("TestIdeals", Reload => true)


p = 7;
m = 1; --m = 1 means working over Z/p^2
e = 3;
R = ZZ[x,y]/(p^(m+1));
f = x^2+y^3;

--p = 7;
--m = 0; --m = 1 means working over Z/p^2
--e = 4;
--R = ZZ[x,y,z]/(p^(m+1));
--f = x*y+z^2;




cartierIdeal = (p, m, e, f, R) ->(
    n := numgens R;
    Rvars := R_*;
    Y := local Y;
    T := ZZ(monoid[(Rvars | toList(Y_1..Y_n)), MonomialOrder=>ProductOrder{n,n},MonomialSize=>64]);
    S:= T/p^(m+1);
    Svars := S_*;
    J := ideal(apply(n,i->Svars#(n+i) - Svars#i^(p^e)))*S;
    h := (substitute(f,S)) % J;
    L:=ideal((coefficients(h,Variables => Rvars))#1);
    L = first entries mingens L;
    subRelations := apply(n,i->Svars#(n+i) => Svars#i);
    L = apply(L, g ->substitute(g,subRelations));
    substitute(ideal L, R)
); -- This function calculates the ideal C^e.f


findJumps = (p,m,e,f,R) ->(
    g := f;
    ideal1 := cartierIdeal(p, m, e, substitute(1,R), R);
    ideal2 := cartierIdeal(p, m, e, g, R);
    jumps := {};
    for s from 0 to p^(e+m)-1 do(
	if ideal1 != ideal2 then(
	    jumps = append(jumps, s);
	    );
    	g = f*g;
    	ideal1 = ideal2;
    	ideal2 = cartierIdeal(p,m,e,g,R);
	print s;
    	);
    jumps
); --- this function is more direct, and perhaps a bit faster, than the creation of ideals list and  computation of jumps done above

cI = new MutableList from toList (p^(e+m)+1:ideal(0_R));
fRaised = new MutableList from toList (p^(e+m)+1:0_R);

-- Precompute all the powers
fRaised#0=1_R;
for i from 0 to p^(e+m) do(
    fRaised#(i+1) = f*fRaised#i;
);


cartierIdealPow = (p, m, e, a, f, R) ->(
    --print(cI#a);
    --print(ideal(0_R));
    --print(cI#a != ideal(0_R));
    if cI#a != ideal(0_R) then return cI#a;
    n := numgens R;
    Rvars := R_*;
    Y := local Y;
    T := ZZ(monoid[(Rvars | toList(Y_1..Y_n)), MonomialOrder=>ProductOrder{n,n},MonomialSize=>64]);
    S:= T/p^(m+1);
    Svars := S_*;
    J := ideal(apply(n,i->Svars#(n+i) - Svars#i^(p^e)))*S;
    h := (substitute(fRaised#a,S)) % J;
    L:=ideal((coefficients(h,Variables => Rvars))#1);
    L = first entries mingens L;
    subRelations := apply(n,i->Svars#(n+i) => Svars#i);
    L = apply(L, g ->substitute(g,subRelations));
    cI#a = substitute(ideal L, R);
    return cI#a;
); -- This function calculates the ideal C^e.f




findJumpsbtw = (a,b,p,m,e,f,R) ->(
    print (a,b);
    if cartierIdealPow(p,m,e,a,f, R) == cartierIdealPow(p,m,e,b,f,R) then return {};
    if a+1 == b then return {a};
    return join( findJumpsbtw(a,(a+b)//2, p, m, e, f, R), findJumpsbtw((a+b)//2, b, p, m, e, f, R) );
);



jumps = findJumpsbtw(1, p^(e+m), p,m,e,f,R);
print(jumps);

file = "calculationsheet.txt";
file << "c" << endl << close;





for g in jumps do(
    print g;
    print adicExpansion(p, g);
    print(trim(cartierIdealPow(p,m,e,g,f,R)));
    print "----";
    );
end

