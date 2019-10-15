# Proof-of-concept implementation of https://github.com/monero-project/research-lab/issues/56
# This version proves knowledge of one or more commitments to zero in a list, with key images and amount commitments

from dumb25519 import hash_to_point, random_scalar, Scalar, hash_to_scalar, G, random_point, Z
import dumb25519
import random

H = hash_to_point('H')

# Proof structure
class Proof:
    def __init__(self):
        self.J = None # key images
        self.A = None
        self.B = None
        self.C = None
        self.D = None
        self.X = None # for signing keys
        self.Y = None # for key images
        self.Z = None # for amounts
        self.f = None
        self.zA = None
        self.zC = None
        self.zR = None # for signing keys and key images
        self.zS = None # for amounts

    def __repr__(self):
        temp = '<TriptychProof> '
        temp += 'J:'+repr(self.J)+'|'
        temp += 'A:'+repr(self.A)+'|'
        temp += 'B:'+repr(self.B)+'|'
        temp += 'C:'+repr(self.C)+'|'
        temp += 'D:'+repr(self.D)+'|'
        temp += 'X:'+repr(self.X)+'|'
        temp += 'Y:'+repr(self.Y)+'|'
        temp += 'Z:'+repr(self.Z)+'|'
        temp += 'f:'+repr(self.f)+'|'
        temp += 'zA:'+repr(self.zA)+'|'
        temp += 'zC:'+repr(self.zC)+'|'
        temp += 'zR:'+repr(self.zR)+'|'
        temp += 'zS:'+repr(self.zS)
        return temp

# Pedersen vector commitment
def com_tensor(v,r):
    C = dumb25519.Z
    for i in range(len(v)):
        for j in range(len(v[0])):
            for k in range(len(v[0][0])):
                C += hash_to_point('Gi',i,j,k)*v[i][j][k]
    C += r*H
    return C

# Decompose a value with given base and size
def decompose(val,base,size,t=int):
    temp_val = val # for reconstruction testing
    r = []
    for i in range(size-1,-1,-1):
        slot = base**i
        r.append(int(val/slot))
        val -= slot*r[-1]
    r = list(reversed(r))

    # Reconstruct to ensure correct decomposition
    if t == int:
        temp = 0
    else:
        temp = Scalar(0)

    for i in range(len(r)):
        temp += r[i]*base**i
    if not temp == temp_val:
        raise ArithmeticError('Decomposition failed!')

    # Return int list or Scalar list
    if t == int:
        return r

    s = []
    for i in r:
        s.append(Scalar(i))
    return s

# Kronecker delta
def delta(x,y):
    if x == y:
        return Scalar(1)
    return Scalar(0)

# Compute a convolution
def convolve(x,y,size=None):
    r = [Scalar(0)]*(len(x)+len(y))
    for i in range(len(x)):
        for j in range(len(y)):
            r[i+j] += x[i]*y[j]

    # Pad with zeros
    if size is not None and size > len(r):
        for i in range(size-len(r)):
            r.append(Scalar(0))

    return r

# Perform a multi-index commitment-to-zero proof
#
# INPUT
#  M: public key list
#  P: input commitment list
#  Q: output commitment list
#  l: list of indices such that each M[l[i]] is a commitment to zero
#  r: list of Pedersen blinders for all M[l[i]]
#  s: list of Pedersen blinders for all P[l[i]]
#  t: list of Pedersen blinders for all Q[l[i]]
#  a: list of Pedersen values for all P[l[i]]
#  b: list of Pedersen values for all Q[l[i]]
#  m: dimension such that len(M) == 2**m
# RETURNS
#  proof structure
def prove(M,P,Q,l,r,s,t,a,b,m):
    n = 2 # binary decomposition

    # Size check
    if not len(M) == n**m or not len(P) == n**m:
        raise IndexError('Bad size decomposition!')
    N = len(M)

    # Reconstruct the known commitments
    if not len(l) == len(r) or not len(l) == len(s):
        raise IndexError('Index/blinder size mismatch!')
    w = len(l)
    for i in range(w):
        if not M[l[i]] == r[i]*G:
            raise ValueError('Bad commitment blinder!')

    # Construct key images
    J = []
    for i in range(w):
        J.append(r[i].invert()*hash_to_point(M[l[i]]))

    # Prepare matrices and corresponding blinders
    rA = random_scalar()
    rB = random_scalar()
    rC = random_scalar()
    rD = random_scalar()

    # Commit to zero-sum blinders
    a = [[[random_scalar() for _ in range(n)] for _ in range(m)] for _ in range(w)]
    for j in range(m):
        for u in range(w):
            a[u][j][0] = Scalar(0)
        for i in range(1,n):
            for u in range(w):
                a[u][j][0] -= a[u][j][i]
    A = com_tensor(a,rA)

    # Commit to decomposition bits
    decomp_l = []
    for u in range(w):
        decomp_l.append(decompose(l[u],n,m))
    sigma = [[[None for _ in range(n)] for _ in range(m)] for _ in range(w)]
    for j in range(m):
        for i in range(n):
            for u in range(w):
                sigma[u][j][i] = delta(decomp_l[u][j],i)
    B = com_tensor(sigma,rB)

    # Commit to a/sigma relationships
    a_sigma = [[[Scalar(0) for _ in range(n)] for _ in range(m)] for _ in range(w)]
    for j in range(m):
        for i in range(n):
            for u in range(w):
                a_sigma[u][j][i] = a[u][j][i]*(Scalar(1) - Scalar(2)*sigma[u][j][i])
    C = com_tensor(a_sigma,rC)
    
    # Commit to squared a-values
    a_sq = [[[Scalar(0) for _ in range(n)] for _ in range(m)] for _ in range(w)]
    for j in range(m):
        for i in range(n):
            for u in range(w):
                a_sq[u][j][i] = -a[u][j][i]**2
    D = com_tensor(a_sq,rD)

    # Compute p coefficients
    p = [[[Scalar(0) for _ in range(m)] for _ in range(N)] for _ in range(w)]
    for k in range(N):
        decomp_k = decompose(k,n,m)
        for u in range(w):
            p[u][k] = [a[u][0][decomp_k[0]],delta(decomp_l[u][0],decomp_k[0])]
        
        for j in range(1,m):
            for u in range(w):
                p[u][k] = convolve(p[u][k],[a[u][j][decomp_k[j]],delta(decomp_l[u][j],decomp_k[j])],m)

        # Combine to single coefficients in p[0]
        for j in range(m):
            for u in range(1,w):
                p[0][k][j] += p[u][k][j]

    # Generate proof values
    X = [dumb25519.Z for _ in range(m)]
    Y = [dumb25519.Z for _ in range(m)]
    Z = [dumb25519.Z for _ in range(m)]
    rho = [[random_scalar() for _ in range(m)] for _ in range(w)]
    for j in range(m):
        for i in range(N):
            X[j] += M[i]*p[0][i][j]
            Y[j] += hash_to_point(M[i])*p[0][i][j]
            Z[j] += P[i]*p[0][i][j]
        for u in range(w):
            X[j] += rho[u][j]*G
            Y[j] += rho[u][j]*J[u]
            Z[j] += rho[u][j]*G

    # Partial proof
    proof = Proof()
    proof.J = J
    proof.A = A
    proof.B = B
    proof.C = C
    proof.D = D
    proof.X = X
    proof.Y = Y
    proof.Z = Z

    # Fiat-Shamir transcript challenge
    x = hash_to_scalar(M,P,Q,J,A,B,C,D,X,Y,Z)

    f = [[[None for _ in range(n)] for _ in range(m)] for _ in range(w)]
    for j in range(m):
        for i in range(1,n):
            for u in range(w):
                f[u][j][i] = sigma[u][j][i]*x + a[u][j][i]

    zA = rB*x + rA
    zC = rC*x + rD
    zR = []
    zS = Scalar(0)
    for u in range(w):
        zR.append(r[u]*x**m)
        zS += s[u]*x**m
    for j in range(m):
        for u in range(w):
            zR[u] -= rho[u][j]*x**j
            zS -= rho[u][j]*x**j
    for i in range(len(t)):
        zS -= t[i]*x**m

    # Assemble proof
    proof.f = f
    proof.zA = zA
    proof.zC = zC
    proof.zR = zR
    proof.zS = zS

    return proof

# Verify a commitment-to-zero proof
#
# INPUT
#  M: public key list
#  P: input commitment list
#  Q: output commitment list
#  proof: proof structure
#  m: dimension such that len(M) == 2**m
# RETURNS
#  True if the proof is valid
def verify(M,P,Q,proof,m):
    n = 2
    N = n**m

    J = proof.J
    A = proof.A
    B = proof.B
    C = proof.C
    D = proof.D
    X = proof.X
    Y = proof.Y
    Z = proof.Z
    f = proof.f
    zA = proof.zA
    zC = proof.zC
    zR = proof.zR
    zS = proof.zS

    # Fiat-Shamir transcript challenge
    x = hash_to_scalar(M,P,Q,J,A,B,C,D,X,Y,Z)

    w = len(J)

    # A/B check
    for j in range(m):
        for u in range(w):
            f[u][j][0] = x
        for i in range(1,n):
            for u in range(w):
                f[u][j][0] -= f[u][j][i]
    if not com_tensor(f,zA) == B*x + A:
        raise ArithmeticError('Failed A/B check!')

    # C/D check
    fx = [[[None for _ in range(n)] for _ in range(m)] for _ in range(w)]
    for j in range(m):
        for i in range(n):
            for u in range(w):
                fx[u][j][i] = f[u][j][i]*(x-f[u][j][i])
    if not com_tensor(fx,zC) == C*x + D:
        raise ArithmeticError('Failed C/D check!')

    # Commitment check
    RX = dumb25519.Z
    RY = dumb25519.Z
    RZ = dumb25519.Z
    for i in range(N):
        t = [Scalar(1) for _ in range(w)]
        decomp_i = decompose(i,n,m)
        for j in range(m):
            for u in range(w):
                t[u] *= f[u][j][decomp_i[j]]
        for u in range(w):
            RX += M[i]*t[u]
            RY += hash_to_point(M[i])*t[u]
            RZ += P[i]*t[u]

    for j in range(m):
        RX -= X[j]*x**j
        RY -= Y[j]*x**j
        RZ -= Z[j]*x**j
    for u in range(w):
        RX -= zR[u]*G
        RY -= zR[u]*J[u]
    T = dumb25519.Z
    for j in range(len(Q)):
        T += Q[j]

    if not RX == dumb25519.Z:
        raise ArithmeticError('Failed signing key check!')
    if not RY == dumb25519.Z:
        raise ArithmeticError('Failed linking check!')
    if not RZ - T*x**m == zS*G:
        raise ArithmeticError('Failed balance check!')

    return True


#
# Basic test
#

# Parameters
m = 3 # ring size is 2**m
spends = 2 # number of spent inputs
outs = 2 # number of generated outputs

# Set up the test
l = random.sample(range(2**m),spends) # spend indices
r = [random_scalar() for _ in range(spends)] # signing keys
s = [random_scalar() for _ in range(spends)] # input commitment blinders
a = [random_scalar() for _ in range(spends)] # input commitment amounts
t = [random_scalar() for _ in range(outs)] # output commitment blinders
b = [random_scalar() for _ in range(outs)] # output commitment amounts

# Ensure output amounts balance input amounts
b[0] = Scalar(0)
for i in range(len(a)):
    b[0] += a[i]
for i in range(1,len(b)):
    b[0] -= b[i]

# Set all keys and commitments
M = [random_point() for _ in range(2**m)] # possible signing keys
P = [random_point() for _ in range(2**m)] # corresponding commitments
Q = []
for i in range(outs): # output commitments
    Q.append(t[i]*G + b[i]*H)
for i in range(spends): # spend commitments
    M[l[i]] = r[i]*G
    P[l[i]] = s[i]*G + a[i]*H

# Run test
print 'Proving...'
proof = prove(M,P,Q,l,r,s,t,a,b,m)
print 'Verifying...'
verify(M,P,Q,proof,m)
