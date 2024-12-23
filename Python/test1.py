from sympy import Symbol, groebner

class Ring:
    def __init__(self, p : int, m : int, generators : list[Symbol]):
        self.p = p
        self.m = m
        self.generators = generators

class Ideal:
    def __init__(self, R : Ring, generators : list[Symbol]):
        self.R = R
        self.generators = generators

    def intersection(I : "Ideal", J : "Ideal"):
        if not isinstance(I, Ideal) or not isinstance(J,Ideal):
            raise TypeError("The arguments must be an Ideal")
        if I.R != J.R:
            raise ValueError("The arguments must be of the same ring")

        G = I.generators + J.generators
        G = groebner(G, *I.R.generators, domain='ZZ' if R.p == 0 else f'GF({R.p**R.m})')
        return Ideal(I.R, list(G))
    
    def F(f : Symbol, e : int):
        monoms = f.args
        for i,m in enumerate(monoms):
            monoms[i] = m**(I.R.p**e)


    def fpower(I : "Ideal", e : int): #compute the frobenius power
        for i,f in enumerate(I.generators):
            I.generators[i] = I.generators[i]**(I.R.p**e)
        return Ideal(I.G)



# Example usage:
x, y = Symbol('x'), Symbol('y')
R = Ring(p=7, generators=[x, y])

I1 = Ideal(R, [x**2, x, 3])
I2 = Ideal(R, [y])

I_inter = Ideal.intersection(I1, I2)

print(I_inter.generators)

print(Ideal.fpower(I1, 2).generators)






#R = QQ.old_poly_ring(x,y)

#I = R.ideal(x**2, y)

#print(R/I)
