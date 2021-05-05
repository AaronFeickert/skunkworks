# Proof-of-concept demonstrating a failure of the hardness assumption of Arcturus: https://eprint.iacr.org/2020/312
# Uses the method described here: https://github.com/UkoeHB/break-dual-target-dl

import arcturus
import unittest
import dumb25519
from dumb25519 import random_scalar, random_point, Scalar
import random

G = dumb25519.G
H = arcturus.H
U = arcturus.U

class TestBreak(unittest.TestCase):
    # Produce a false proof by breaking the soundness hardness assumption
    def test_break(self):
        # Proof parameters
        lg_N = 3
        N = 2**lg_N
        spends = 3
        outs = 2

        # The adversary's inputs
        l = random.sample(range(N),spends) # spend indices
        r = [random_scalar() for _ in range(spends)] # signing keys
        s = [random_scalar() for _ in range(spends)] # input commitment blinders
        a = [random_scalar() for _ in range(spends)] # input commitment amounts
        t = [random_scalar() for _ in range(outs)] # output commitment blinders
        b = [random_scalar() for _ in range(outs)] # output commitment amounts

        # Balance amounts
        b[0] = Scalar(0)
        for i in range(len(a)):
            b[0] += a[i]
        for i in range(1,len(b)):
            b[0] -= b[i]

        # Set keys and commitments
        M = [random_point() for _ in range(N)]
        P = [random_point() for _ in range(N)]
        Q = []
        for i in range(outs): # output commitments
            Q.append(t[i]*G + b[i]*H)
        for i in range(spends): # spend commitments
            M[l[i]] = r[i]*G
            P[l[i]] = s[i]*G + a[i]*H

        # Generate invalid linking tags
        evil_j = [random_scalar() for _ in range(spends)]
        J = [evil_j[i]*U for i in range(spends)]

        # Initial proof
        proof,state = arcturus.prove_initial(M,P,Q,J,l,lg_N)

        # Modify signing data to account for invalid linking tags
        mu = state.mu
        evil_r = [random_scalar() for _ in range(spends)] # for signing keys

        lambda_r = mu**l[0]*r[0] + mu**l[1]*r[1]
        for i in range(2,spends):
            lambda_r += mu**l[i]*(r[i] - evil_r[i])

        lambda_j = mu**l[0] + mu**l[1]
        for i in range(2,spends):
            lambda_j += mu**l[i]*(Scalar(1) - evil_r[i]*evil_j[i])
        
        evil_r[0] = -(lambda_j - evil_j[1]*lambda_r)*(mu**l[0]*(evil_j[1] - evil_j[0])).invert()
        evil_r[1] = (lambda_j - evil_j[0]*lambda_r)*(mu**l[1]*(evil_j[1] - evil_j[0])).invert()

        # Final proof
        proof = arcturus.prove_final(l,evil_r,s,t,proof,state)

        # This verification will succeed, but should fail if soundness holds!
        arcturus.verify(M,P,Q,proof,lg_N)

unittest.TextTestRunner(verbosity=2,failfast=True).run(unittest.TestLoader().loadTestsFromTestCase(TestBreak))
