loadPackage("Dmodules", Reload=>true);

R = QQ[x,y,z, dx, dy, dz, WeylAlgebra => {x=>dx, y=>dy, z=>dz}];
f = x*y + z^2;
P = globalBoperator(f);
b = globalBFunction(f);
b = factor(b)

print(b);
print(P);
