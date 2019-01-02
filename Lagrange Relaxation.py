# conda install -c r rpy2

# load the sympy package

from sympy import *
from sympy.simplify import simplify as _simplify, signsimp, nsimplify
import sympy as sp

# create the following functions

# the following function will compute the hessian matrix of any differentiable function

def hessian(f, varlist):

    if isinstance(varlist, (list, tuple)):
        m = len(varlist)
    elif isinstance(varlist, Matrix):
        m = varlist.cols
        assert varlist.rows == 1
    else:
        raise ValueError("Improper variable list in hessian function")
    assert m > 0
    try:
        f.diff(varlist[0]) 
    except AttributeError:
        raise ValueError("Function %d is not differentiable" % i)
    out = zeros(m)
    for i in range(m):
        for j in range(i,m):
            out[i,j] = f.diff(varlist[i]).diff(varlist[j])
    for i in range(m):
        for j in range(i):
            out[i,j] = out[j,i]
    return out

# the following function runs the berkowitz algorithm
    # this algorithm is useful for computing derminants, principal minors, eigenvalues, and the characteristic polynomial of an N x N matrix

def berkowitz(self):

        from sympy.matrices import zeros
        berk = ((1,),)
        if not self:
            return berk

        if not self.is_square:
            raise NonSquareMatrixError()

        A, N = self, self.rows
        transforms = [0]*(N - 1)

        for n in range(N, 1, -1):
            T, k = zeros(n + 1, n), n - 1

            R, C = -A[k, :k], A[:k, k]
            A, a = A[:k, :k], -A[k, k]

            items = [C]

            for i in range(0, n - 2):
                items.append(A*items[i])

            for i, B in enumerate(items):
                items[i] = (R*B)[0, 0]

            items = [S.One, a] + items

            for i in range(n):
                T[i:, i] = items[:n - i + 1]

            transforms[k - 1] = T

        polys = [self._new([S.One, -A[0, 0]])]

        for i, T in enumerate(transforms):
            polys.append(T*polys[i])

        return berk + tuple(map(tuple, polys))

# this function extracts the determinant from the berkowitz algorithm

def berkowitz_det(self):

	if not self.is_square:
		raise NonSquareMatrixError()
	if not self:
		return S.One
	poly = self.berkowitz()[-1]
	sign = (-1)**(len(poly) - 1)
	return sign*poly[-1]

# this function extracts the characteristic polynomial from the berkowitz algorithm

def berkowitz_charpoly(self, x=Dummy('lambda'), simplify=_simplify):

	return PurePoly(list(map(simplify, self.berkowitz()[-1])), x)

# this function extracts the eigenvalues from the berkowitz algorithm

def berkowitz_eigenvals(self, **flags):

	return roots(self.berkowitz_charpoly(Dummy('x')), **flags)	

# this function extracts the principal minors from the berkowitz algorithm

def berkowitz_minors(self):

    sign, minors = S.One, []

    for poly in self.berkowitz():
        minors.append(sign*poly[-1])
        sign = -sign

    return tuple(minors)

# how to do differentiation

    # declare variables
	
x = Symbol('x')
y = Symbol('y')

    # single derivatives
	
diff(sin(y), 
     y)

    # partial derivatives
	
diff(x * sin(y), 
     x)

diff(x * sin(y), 
     y)

    # second order derivatives
	
diff(x * sin(y), 
     y, 2)

diff(x * sin(y), 
     x, 2)

    # second order partial derivatives
        # y then x
		
diff(x * sin(y), 
     y, x)
        
        # x then y
		
diff(x * sin(y), 
     x, y)

# how to compute the hessian

fun = x * sin(y)
mylist = [x, y]

mat = hessian(fun, mylist)

mat

# how to compute the derminant, principal minors, eigenvalues, and characteristic polynomial

berkowitz_det(mat)

berkowitz_minors(mat)

berkowitz_eigenvals(mat)

berkowitz_charpoly(mat)

# now lets compute the hessian for the obj function from NLP.mod

    # define all variables

x11 = Symbol('x11', real = True, nonnegative = True)
x21 = Symbol('x21', real = True, nonnegative = True)
x31 = Symbol('x31', real = True, nonnegative = True)
x41 = Symbol('x41', real = True, nonnegative = True)
x12 = Symbol('x12', real = True, nonnegative = True)
x22 = Symbol('x22', real = True, nonnegative = True)
x32 = Symbol('x32', real = True, nonnegative = True)
x42 = Symbol('x42', real = True, nonnegative = True)
x13 = Symbol('x13', real = True, nonnegative = True)
x23 = Symbol('x23', real = True, nonnegative = True)
x33 = Symbol('x33', real = True, nonnegative = True)
x43 = Symbol('x43', real = True, nonnegative = True)
p1 = Symbol('p1', real = True, nonnegative = True)
p2 = Symbol('p2', real = True, nonnegative = True)
p3 = Symbol('p3', real = True, nonnegative = True)
p4 = Symbol('p4', real = True, nonnegative = True)

    # define objective function, and the variables in a list
    
obj = (12.4*x11) + ((1 - (0.00325884491899048*(x11 / 400)))*23400*x11) + (28.7*x21) + ((1 - (0.00557401772798191*(x21 / 300)))*15100*x21) + (21.4*x31) + ((1 - (0.00233690127478853*(x31 / 540)))*11600*x31) + (22*x41) + ((1 - (0.00412312138699061*(x41 / 150)))*20500*x41) + (22.9*x12) + ((1 - (0.00325884491899048*(x12 / 400)))*23400*x12) + (27.5*x22) + ((1 - (0.00557401772798191*(x22 / 300)))*15100*x22) + (27.4*x32) + ((1 - (0.00233690127478853*(x32 / 540)))*11600*x32) + (28.4*x42) + ((1 - (0.00412312138699061*(x42 / 150)))*20500*x42) + (10.1*x13) + ((1 - (0.00325884491899048*(x13 / 400)))*23400*x13) + (21.6*x23) + ((1 - (0.00557401772798191*(x23 / 300)))*15100*x23) + (30.5*x33) + ((1 - (0.00233690127478853*(x33 / 540)))*11600*x33) + (33.7*x43) + ((1 - (0.00412312138699061*(x43 / 150)))*20500*x43)

varlist = [x11, x21, x31, x41, x12, x22, x32, x42, x13, x23, x33, x43, p1, p2, p3, p4]

    # compute the hessian
    
obj_H = hessian(obj,
                varlist)

obj_H

# export obj_H
    # lets figure out where our current work directory is
    # then we will know how to specify the file path to our obj.txt file

import os
os.getcwd()

    # lets compute the number of rows of the matrix
    # this is so we know how we will need to loop when writing to a text file

len(obj_H[:,1]) # total number of rows

    # open the location of obj.txt file
    # 'w' specifies that we are writing to a file, which will overwrite anything currently in that file
  
objtxt = open('C:\\Users\\Nick\\Desktop\\obj.txt', 'w')

    # loop through every row in obj_H
    # convert every row into a list, then into a string, then add a '\n' to the end so the text file starts a new row before writing the next row
    
for row in range(0, len(obj_H[:,1])):
    objtxt.write(str(str(list(obj_H[[row],:])) + '\n'))

    # close off the location to the obj.txt file
    
objtxt.close()

    # lets compute the principal minors of the hessian

obj_minors = berkowitz_minors(obj_H)

obj_minors

    # lets append these minors to the obj.txt file
    # 'a' specifies that we are appending to a file, which will add to anything currently in that file
    # we convert obj_minors to a string, then add '\n' to the start so the text file knows to start a new line before writing obj_minors 
    
objtxt = open('C:\\Users\\Nick\\Desktop\\obj.txt', 'a')
objtxt.write(str('\n' + str(obj_minors)))
objtxt.close()

# all constriants are equalities and linear therefore they are all convex and differentiable
    # so we don't need to compute the hessians principal minors of each constriant
    # now that we know the objective function and feasible region are both convex, and this is a minimization problem, we know that any stationary point is optimal
    # we also now know that when we perform largrange relaxation, that the resulting function is convex becuase it is a combination of convex functions
    
# lets solve for the values of every variable by using Lagrange relaxation

    # define penalty variables

lam11 = Symbol('lam11', real = True)
lam12 = Symbol('lam12', real = True)
lam13 = Symbol('lam13', real = True)
lam21 = Symbol('lam21', real = True)
lam22 = Symbol('lam22', real = True)
lam23 = Symbol('lam23', real = True)
lam24 = Symbol('lam24', real = True)

    # define relaxed obj

obj_relax = (12.4*x11) + ((1 - (0.00325884491899048*(x11 / 400)))*23400*x11) + (28.7*x21) + ((1 - (0.00557401772798191*(x21 / 300)))*15100*x21) + (21.4*x31) + ((1 - (0.00233690127478853*(x31 / 540)))*11600*x31) + (22*x41) + ((1 - (0.00412312138699061*(x41 / 150)))*20500*x41) + (22.9*x12) + ((1 - (0.00325884491899048*(x12 / 400)))*23400*x12) + (27.5*x22) + ((1 - (0.00557401772798191*(x22 / 300)))*15100*x22) + (27.4*x32) + ((1 - (0.00233690127478853*(x32 / 540)))*11600*x32) + (28.4*x42) + ((1 - (0.00412312138699061*(x42 / 150)))*20500*x42) + (10.1*x13) + ((1 - (0.00325884491899048*(x13 / 400)))*23400*x13) + (21.6*x23) + ((1 - (0.00557401772798191*(x23 / 300)))*15100*x23) + (30.5*x33) + ((1 - (0.00233690127478853*(x33 / 540)))*11600*x33) + (33.7*x43) + ((1 - (0.00412312138699061*(x43 / 150)))*20500*x43) + (lam11 * (x11 + x21 + x31 + x41 - 1500)) + (lam12 * (x12 + x22 + x32 + x42 - 2750)) + (lam13 * (x13 + x23 + x33 + x43 - 3000)) + (lam21 * (x11 + x12 + x13 + p1 - 2000)) + (lam22 * (x21 + x22 + x23 + p2 - 1800)) + (lam23 * (x31 + x32 + x33 + p3 - 2700)) + (lam24 * (x41 + x42 + x43 + p4 - 1000))

    # redefine the variable list

varlist = [x11, x21, x31, x41, x12, x22, x32, x42, x13, x23, x33, x43, p1, p2, p3, p4, lam11, lam12, lam13, lam21, lam22, lam23, lam24]

    # compute the gradient
    
obj_relax_grad = [sp.Eq(diff(obj_relax, i)) for i in varlist]

    # verify that the following objects are lists
                  
type(obj_relax_grad) is list

type(varlist) is list

    # try to solve the system of linear equations
    
linsolve(obj_relax_grad, varlist)

    # lets try the solve command instead

sp.solve(obj_relax_grad, varlist)

# the linsolve doesn't respect the nonnegativity defined for the symbols 
# lets solve by creating an linear program, using the obj_grad equations as contraints
# write obj_grad to a text file

gradtxt = open('C:\\Users\\Nick\\Desktop\\grad.txt', 'w')

for eqn in range(0, len(obj_relax_grad)):
    gradtxt.write(str(str(obj_relax_grad[eqn]) + '\n'))
    
gradtxt.close()



