from common import *
from dumb25519 import *

# Internal proof state
class State:
    def __init__(self):
        self.m = None
        self.n = None
        self.sigma = None
        self.a = None
        self.rA = None
        self.rB = None
        self.rC = None
        self.rD = None
        self.v = None
        self.r = None
        self.rho = None
        self.tau = None
        self.gammas = None

        self.x = None # Fiat-Shamir challenge

# Recursion data
class RecursionData:
    def __init__(self,n,m,f,scalars):
        self.begin = 0
        self.end = n**m
        self.n = n
        self.f = f
        self.scalars = scalars

# Proof structure
class Proof:
    def __init__(self):
        self.A = None
        self.B = None
        self.C = None
        self.D = None
        self.G = None
        self.Q = None
        self.f = None
        self.zA = None
        self.zC = None
        self.zV = None
        self.zR = None

    def __repr__(self):
        temp = '<GrothProof> '
        temp += 'A:'+repr(self.A)+'|'
        temp += 'B:'+repr(self.B)+'|'
        temp += 'C:'+repr(self.C)+'|'
        temp += 'D:'+repr(self.D)+'|'
        temp += 'G:'+repr(self.G)+'|'
        temp += 'Q:'+repr(self.Q)+'|'
        temp += 'f:'+repr(self.f)+'|'
        temp += 'zA:'+repr(self.zA)+'|'
        temp += 'zC:'+repr(self.zC)+'|'
        temp += 'zV:'+repr(self.zV)+'|'
        temp += 'zR:'+repr(self.zR)
        return temp

# Compute an aggregate Fiat-Shamir challenge
def challenge(M,proofs):
    A = [proof.A for proof in proofs]
    B = [proof.B for proof in proofs]
    C = [proof.C for proof in proofs]
    D = [proof.D for proof in proofs]
    G = [proof.G for proof in proofs]
    Q = [proof.Q for proof in proofs]
    return hash_to_scalar(M,A,B,C,D,G,Q)

# Double-blinded Pedersen matrix commitment
def com_matrix(v,r):
    C = r*H
    for j in range(len(v)):
        for i in range(len(v[0])):
            C += hash_to_point('Gi',j,i)*v[j][i]
    return C

# Double-blinded Pedersen commitment
def comm(s,v,r):
    return G*s + H1*v + H2*r

# Recursively build commitment scalars
def recurse(t_,j,data):
    j -= 1
    if j == -1:
        if data.begin < data.end:
            data.scalars[data.begin] = t_
            data.begin += 1
        return
    
    for i in range(data.n):
        t = data.f[j][i]
        t *= t_
        recurse(t,j,data)

# Decompose a value with given base and size
def decompose(val,base,size):
    r = []
    for i in range(size-1,-1,-1):
        slot = base**i
        r.append(val//slot)
        val -= slot*r[-1]
    r = list(reversed(r))

    return r

# Kronecker delta
def delta(x,y):
    if x == y:
        return Scalar(1)
    return Scalar(0)

# Compute a convolution with a degree-one polynomial
def convolve(x,y):
    if not len(y) == 2:
        raise ValueError('Convolution requires a degree-one polynomial!')

    r = [Scalar(0)]*(len(x)+1)
    for i in range(len(x)):
        for j in range(len(y)):
            r[i+j] += x[i]*y[j]

    return r

# Perform a commitment-to-zero proof
#
# INPUT
#  M: list of double-blinded Pedersen commitments such that len(M) == n**m
#  l: index such that M[l] is a commitment to zero
#  v: first Pedersen blinder for M[l]
#  r: second Pedersen blinder for M[l]
#  n,m: dimensions such that len(M) == n**m
# RETURNS
#  proof structure
#  internal state
def prove_initial(M,l,v,r,n,m):
    # Size check
    if not n > 1:
        raise ValueError('Must have nontrivial decomposition base!')
    if not m > 1:
        raise ValueError('Must have m > 1 as decomposition exponent!')
    if not len(M) == n**m:
        raise IndexError('Bad size decomposition!')
    N = len(M)

    # Reconstruct the known commitment
    if not comm(Scalar(0),v,r) == M[l]:
        raise ValueError('Bad known commitment!')

    rA = random_scalar()
    rB = random_scalar()
    rC = random_scalar()
    rD = random_scalar()

    # Commit to zero-sum blinders
    a = [[random_scalar() for _ in range(n)] for _ in range(m)]
    for j in range(m):
        a[j][0] = Scalar(0)
        for i in range(1,n):
            a[j][0] -= a[j][i]
    A = com_matrix(a,rA)

    # Commit to decomposition bits
    decomp_l = decompose(l,n,m)
    sigma = [[None for _ in range(n)] for _ in range(m)]
    for j in range(m):
        for i in range(n):
            sigma[j][i] = delta(decomp_l[j],i)
    B = com_matrix(sigma,rB)

    # Commit to a/sigma relationships
    a_sigma = [[Scalar(0) for _ in range(n)] for _ in range(m)]
    for j in range(m):
        for i in range(n):
            a_sigma[j][i] = a[j][i]*(Scalar(1) - Scalar(2)*sigma[j][i])
    C = com_matrix(a_sigma,rC)
    
    # Commit to squared a-values
    a_sq = [[Scalar(0) for _ in range(n)] for _ in range(m)]
    for j in range(m):
        for i in range(n):
            a_sq[j][i] = -a[j][i]*a[j][i]
    D = com_matrix(a_sq,rD)

    # Compute p coefficients
    p = [[] for _ in range(N)]
    for k in range(N):
        decomp_k = decompose(k,n,m)
        p[k] = [a[0][decomp_k[0]],delta(decomp_l[0],decomp_k[0])]
        
        for j in range(1,m):
            p[k] = convolve(p[k],[a[j][decomp_k[j]],delta(decomp_l[j],decomp_k[j])])

    # Generate proof values
    G = [Z for _ in range(m)]
    Q = [Z for _ in range(m)]
    rho = [None for _ in range(m)]
    tau = [None for _ in range(m)]
    gammas = [None for _ in range(m)]
    for j in range(m):
        rho[j] = random_scalar()
        tau[j] = random_scalar()
        gamma = random_scalar()
        gammas[j] = gamma
        for i in range(N):
            G[j] += M[i]*p[i][j]
        G[j] -= H2*gamma
        Q[j] = comm(Scalar(0),rho[j],tau[j]) + H2*gamma

    # Assemble state
    state = State()
    state.m = m
    state.n = n
    state.sigma = sigma
    state.a = a
    state.rA = rA
    state.rB = rB
    state.rC = rC
    state.rD = rD
    state.v = v
    state.r = r
    state.rho = rho
    state.tau = tau
    state.gammas = gammas

    # Partial proof
    proof = Proof()
    proof.A = A
    proof.B = B
    proof.C = C
    proof.D = D
    proof.G = G
    proof.Q = Q

    return proof,state

# Complete a partial proof
def prove_final(proof,state):
    x = state.x # aggregate Fiat-Shamir challenge

    # Recover state
    m = state.m
    n = state.n
    sigma = state.sigma
    a = state.a
    rA = state.rA
    rB = state.rB
    rC = state.rC
    rD = state.rD
    v = state.v
    r = state.r
    rho = state.rho
    tau = state.tau
    gammas = state.gammas

    f = [[None for _ in range(n-1)] for _ in range(m)]
    for j in range(m):
        for i in range(1,n):
            f[j][i-1] = sigma[j][i]*x + a[j][i]

    zA = rB*x + rA
    zC = rC*x + rD
    zV = v*x**m
    zR = r*x**m
    for j in range(m):
        zV -= rho[j]*x**j
        zR -= tau[j]*x**j

    # Assemble proof
    proof.f = f
    proof.zA = zA
    proof.zC = zC
    proof.zV = zV
    proof.zR = zR

    return proof,gammas

# Verify a commitment-to-zero proof
#
# INPUT
#  M: list of double-blinded Pedersen commitments such that len(M) == n**m
#  proof: proof structure
#  n,m: dimensions such that len(M) == n**m
#  x: aggregate Fiat-Shamir challenge (in practice this is computed by the verifier)
# RETURNS
#  True if the proof is valid
def verify(M,proof,n,m,x):
    A = proof.A
    B = proof.B
    C = proof.C
    D = proof.D
    G = proof.G
    Q = proof.Q
    f = [[None for _ in range(n)] for _ in range(m)]
    zA = proof.zA
    zC = proof.zC
    zV = proof.zV
    zR = proof.zR

    N = n**m

    for j in range(m):
        f[j][0] = x
        for i in range(1,n):
            f[j][i] = proof.f[j][i-1]
            f[j][0] -= f[j][i]

    # A/B check
    if not com_matrix(f,zA) == B*x + A:
        raise ArithmeticError('Failed A/B check!')

    # C/D check
    fx = [[None for _ in range(n)] for _ in range(m)]
    for j in range(m):
        for i in range(n):
            fx[j][i] = f[j][i]*(x-f[j][i])
    if not com_matrix(fx,zC) == C*x + D:
        raise ArithmeticError('Failed C/D check!')

    # Commitment check
    scalars = ScalarVector([])
    points = PointVector([])
    s = ScalarVector([Scalar(0) for i in range(N)])
    recurse(Scalar(1),m,RecursionData(n,m,f,s))
    for i in range(N):
        scalars.append(s[i])
        points.append(M[i])
    for j in range(m):
        scalars.append(-x**j)
        points.append(G[j])
        scalars.append(-x**j)
        points.append(Q[j])
    if not multiexp(scalars,points) == comm(Scalar(0),zV,zR):
        raise ArithmeticError('Failed commitment check!')

    return True
