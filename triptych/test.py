import triptych
import unittest
import dumb25519
from dumb25519 import hash_to_scalar, random_scalar, random_point, Scalar
import random

G = dumb25519.G
H = dumb25519.hash_to_point('H')

class TestValidProofs(unittest.TestCase):
    def test_valid_prove_verify(self):
        print('')

        for m in [2,3]: # ring size 4,8
            for w in [1,2]: # number of proofs in batch
                print('Test parameter (m,w):',m,w)
                l = random.sample(range(2**m),w) # spend indices
                r = [random_scalar() for _ in range(w)] # signing key
                s = [random_scalar() for _ in range(w)] # commitment key

                # Data to hide
                seeds = [random_scalar() for _ in range(w)] # used to derive the blinding factors for the verifier; in this implementation, must be unique across proofs
                offsets = [random_scalar() for _ in range(w)] # offset to hide the embedded values from the verifier
                blinders = [[hash_to_scalar('first blinder',seeds[u]),hash_to_scalar('second blinder',seeds[u])] for u in range(w)]
                aux = [random_scalar() for _ in range(w)] # data to embed into proofs for later recovery

                # Set keys and run proofs
                M = [random_point() for _ in range(2**m)] # possible signing keys
                P = [random_point() for _ in range(2**m)] # corresponding commitments
                proofs = []
                for u in range(w):
                    M[l[u]] = r[u]*G
                    P[l[u]] = s[u]*G
                for u in range(w):
                    proof = triptych.prove(M,P,l[u],r[u],s[u],m,blinders[u],aux[u])
                    proof.embeds = [blinders[u][0],blinders[u][1] + offsets[u]]
                    proofs.append(proof)

                # Verify all proofs in batch
                offset_aux = triptych.verify(M,P,proofs,m)
                for u in range(w):
                    self.assertEqual(aux[u],offset_aux[u] + offsets[u])

unittest.TextTestRunner(verbosity=2,failfast=True).run(unittest.TestLoader().loadTestsFromTestCase(TestValidProofs))
