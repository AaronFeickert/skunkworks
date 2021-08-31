# Coin structure

from dumb25519 import Point, Scalar, PointVector, ScalarVector, random_scalar, hash_to_scalar
import address
import bpplus

class CoinParameters:
	def __init__(self,G,F,H,N):
		if not isinstance(G,Point):
			raise TypeError('Bad type for parameter G!')
		if not isinstance(F,Point):
			raise TypeError('Bad type for parameter F!')
		if not isinstance(H,Point):
			raise TypeError('Bad type for parameter H!')
		if not isinstance(N,int) or N < 1:
			raise ValueError('Bad type or value for parameter N!')
		
		self.G = G
		self.F = F
		self.H = H
		self.N = N

class Coin:
	def __init__(self,params,public,value,memo,is_output=False):
		if not isinstance(params,CoinParameters):
			raise TypeError('Bad type for parameters!')
		if not isinstance(public,address.PublicAddress):
			raise TypeError('Bad type for coin address!')
		if not isinstance(value,int) or value < 0 or value >= 2**params.N:
			raise ValueError('Bad type or value for coin value!')
		if not isinstance(memo,str):
			raise TypeError('Bad type for coin memo!')
		if not isinstance(is_output,bool):
			raise TypeError('Bad type for coin output flag!')

		# Recovery key
		k = random_scalar()
		self.K = k*params.G
		K_der = k*public.Q1

		# Serial number commitment
		self.S = hash_to_scalar('ser',K_der)*params.G + public.Q2

		# Value commitment
		self.C = Scalar(value)*params.G + hash_to_scalar('val',K_der)*params.F
		self.range = bpplus.prove(
			bpplus.RangeStatement(bpplus.RangeParameters(params.F,params.G,params.N),PointVector([self.C])),
			bpplus.RangeWitness(ScalarVector([Scalar(value)]),ScalarVector([hash_to_scalar('val',K_der)]))
		)

		# Value and memo
		# NOTE: These are not yet encrypted!
		self.value = Scalar(value)
		self.memo = memo

		# Data used for output only
		self.is_output = False
		if is_output:
			self.is_output = True
			self.k = k
			self.Q1 = public.Q1

		self.recovered = False
	
	def identify(self,params,public,incoming):
		if not isinstance(params,CoinParameters):
			raise TypeError('Bad type for parameters!')
		if not isinstance(public,address.PublicAddress):
			raise TypeError('Bad type for public address!')
		if not isinstance(incoming,address.IncomingViewKey):
			raise TypeError('Bad type for incoming view key!')
	
		K_der = incoming.s1*self.K
		
		# Test for ownership
		if not self.S == hash_to_scalar('ser',K_der)*params.G + public.Q2:
			raise ArithmeticError('Coin does not belong to this public address!')
		
		# Test for value commitment
		if not self.C == self.value*params.G + hash_to_scalar('val',K_der)*params.F:
			raise ArithmeticError('Bad coin value commitment!')
		
		# Test range proof
		bpplus.verify(
			[bpplus.RangeStatement(bpplus.RangeParameters(params.F,params.G,params.N),PointVector([self.C]))],
			[self.range]
		)
		
	def recover(self,params,public,full):
		if not isinstance(params,CoinParameters):
			raise TypeError('Bad type for coin parameters!')
		if not isinstance(public,address.PublicAddress):
			raise TypeError('Bad type for public address!')
		if not isinstance(full,address.FullViewKey):
			raise TypeError('Bad type for full view key!')
		
		K_der = full.s1*self.K
		
		# Test for ownership
		if not self.S == hash_to_scalar('ser',K_der)*params.G + public.Q2:
			raise ArithmeticError('Coin does not belong to this public address!')
		
		# Test for value commitment
		if not self.C == self.value*params.G + hash_to_scalar('val',K_der)*params.F:
			raise ArithmeticError('Bad coin value commitment!')
		
		# Test range proof
		bpplus.verify(
			[bpplus.RangeStatement(bpplus.RangeParameters(params.F,params.G,params.N),PointVector([self.C]))],
			[self.range]
		)
		
		# Recover serial number and generate tag
		self.s = hash_to_scalar('ser',K_der) + full.s2
		self.T = self.s.invert()*params.H

		self.recovered = True
