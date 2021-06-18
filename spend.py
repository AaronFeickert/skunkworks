# Transaction size comparison (spend)
#
# This script compares spend transaction sizes between Lelantus and Triptych

import argparse
import math
import tabulate

# Size of group/scalar elements, in bytes
SIZE = 32

# Parse arguments
parser = argparse.ArgumentParser(description='Compare spend transaction sizes between Lelantus and Triptych')
parser.add_argument('-n',type=int,help='input set size base',required=True)
parser.add_argument('-m',type=int,help='input set size exponent',required=True)
parser.add_argument('-w',type=int,help='number of spends',required=True)
parser.add_argument('-T',type=int,help='number of outputs',required=True)
args = parser.parse_args()

if not (args.n > 1 and args.m > 1 and args.w > 0 and args.T > 0):
	raise argparse.ArgumentError('Bad input argument!')

n = args.n
m = args.m
w = args.w
T = args.T

headers = ['Component', 'Elements', 'Size (B)']

# Lelantus
#
# A spend transaction has the following components, with number of group/scalar elements
# (*) Note: fee is represented by an 8-byte integer
#
# Addresses:				5*T
# Coin commitments:			T
# Encrypted values/keys:	2*T
# Range proof:				2*lg(64*T) + 10
# Serials:					w
# One-of-many proofs:		[m*(n-1) + 2*m + 8]*w
# Balance proof:			T + 2
# Signatures:				2*w
# Fee (*):					{8 bytes}
lelantus = []
lelantus.append(['Addresses', 5*T])
lelantus.append(['Coin commitments', T])
lelantus.append(['Encrypted values/keys', 2*T])
lelantus.append(['Range proof', 2*math.ceil(math.log2(64*T)) + 10])
lelantus.append(['Serials', w])
lelantus.append(['One-of-many proofs', (m*(n - 1) + 2*m + 8)*w])
lelantus.append(['Balance proof', T + 2])
lelantus.append(['Signatures', 2*w])

for i in range(len(lelantus)):
	lelantus[i].append(lelantus[i][1]*SIZE)
lelantus.append(['Fee','',8])

print('Data for Lelantus spend transaction')
print('n = {}, m = {}, w = {}, T = {}'.format(n,m,w,T))
print()
print(tabulate.tabulate(lelantus, headers=headers))
total = sum([row[2] for row in lelantus])
print()
print('Total size is {} kB'.format(round(total/1024,2)))

print()
print()

# Triptych
#
# A spend transaction has the following components, with number of group/scalar elements
# (*) Note: fee is represented by an 8-byte integer
#
# Output keys:				T
# Amount commitments:		T
# Amount offsets:			w
# Range proof:				2*lg(64*T) + 9
# Encrypted values/keys:	2*T
# One-of-many proofs:		[m*(n-1) + 2*m + 8]*w
# Key images:				w
# Diffie-Hellman key:		1
# Fee (*):					{8 bytes}
triptych = []
triptych.append(['Output keys', T])
triptych.append(['Amount commitments', T])
triptych.append(['Amount offsets', w])
triptych.append(['Range proof', 2*math.ceil(math.log2(64*T)) +9])
triptych.append(['Encrypted values/keys', 2*T])
triptych.append(['One-of-many proofs', (m*(n - 1) + 2*m + 8)*w])
triptych.append(['Key images', w])
triptych.append(['Diffie-Hellman key', 1])

for i in range(len(triptych)):
	triptych[i].append(triptych[i][1]*SIZE)
triptych.append(['Fee','',8])

print('Data for Triptych spend transaction')
print('n = {}, m = {}, w = {}, T = {}'.format(n,m,w,T))
print()
print(tabulate.tabulate(triptych, headers=headers))
total = sum([row[2] for row in triptych])
print()
print('Total size is {} kB'.format(round(total/1024,2)))