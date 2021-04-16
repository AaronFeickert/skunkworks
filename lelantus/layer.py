# Two-layer proving system for commitment to zero

from common import *
from dumb25519 import *
import groth
from random import randint

# Total of N = M*T = (n_M**m_M)*(n_T**m_T) commitments
n_M = 2
m_M = 2
n_T = 2
m_T = 2
M = n_M**m_M
T = n_T**m_T
N = M*T

# Generate the commitments
C = PointVector([]) # commitments
L = randint(0,N-1) # secret index
v = random_scalar() # value
r = random_scalar() # mask
for i in range(N):
    C.append(random_point())
C[L] = groth.comm(Scalar(0),v,r)

# Identify the subset major and minor indices for L
k = L // M
l = L % M
if not L == k*M + l:
    raise ArithmeticError('Bad secret index!')

# Prepare offset commitments
masks_v = ScalarVector([])
masks_r = ScalarVector([])
d = PointVector([])
for i in range(M):
    masks_v.append(random_scalar())
    masks_r.append(random_scalar())
    d.append(C[k*M+i] + H1*masks_v[-1] + H2*masks_r[-1])

# Layer-1 proof
proof1,state = groth.prove_initial(d,l,v+masks_v[l],r+masks_r[l],n_M,m_M)
state.x = groth.challenge(d,[proof1])
proof1 = groth.prove_final(proof1,state)[0] # ignore gamma return value for this test

# Fiat-Shamir challenges
prefix = hash_to_scalar(C,d,proof1)
x = ScalarVector([])
for i in range(M):
    x.append(hash_to_scalar(prefix,i))

# Subset digests
D = PointVector([])
for i in range(T):
    D.append(multiexp(x,C[i*M:(i+1)*M]))
digest_d = multiexp(x,d)

# Layer-2 proof
mask_v = Scalar(0)
mask_r = Scalar(0)
for i in range(M):
    mask_v += x[i]*masks_v[i]
    mask_r += x[i]*masks_r[i]
proof2,state = groth.prove_initial([digest_d - D[i] for i in range(T)],k,mask_v,mask_r,n_T,m_T)
state.x = groth.challenge([digest_d - D[i] for i in range(T)],[proof2])
proof2 = groth.prove_final(proof2,state)[0] # ignore gamma return value for this test

# Verify the layer-1 proof
groth.verify(d,proof1,n_M,m_M,groth.challenge(d,[proof1]))

# Verify the layer-2 proof
groth.verify([digest_d - D[i] for i in range(T)],proof2,n_T,m_T,groth.challenge([digest_d - D[i] for i in range(T)],[proof2]))
