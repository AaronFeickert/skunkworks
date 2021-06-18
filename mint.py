# Transaction size comparison (mint)
#
# This script compares mint transaction sizes between Lelantus and Triptych

import argparse
import math
import tabulate

# Size of group/scalar elements, in bytes
SIZE = 32

# Parse arguments
parser = argparse.ArgumentParser(description='Compare mint transaction sizes between Lelantus and Triptych')
parser.add_argument('-T',type=int,help='number of outputs',required=True)
args = parser.parse_args()

if not args.T > 0:
	raise argparse.ArgumentError('Bad input argument!')

T = args.T

headers = ['Component', 'Elements', 'Size (B)']

# Lelantus
#
# A mint transaction has the following components, with number of group/scalar elements
# (*) Note: fee and values are represented by 8-byte integers
#
# Addresses:				5*T
# Coin commitments:			T
# Encrypted keys:			T
# Values (*):				{8 bytes}*T
# Signatures:				2*T
# Fee (*):					{8 bytes}
lelantus = []
lelantus.append(['Addresses', 5*T])
lelantus.append(['Coin commitments', T])
lelantus.append(['Encrypted keys', T])
lelantus.append(['Signatures', 2*T])

for i in range(len(lelantus)):
	lelantus[i].append(lelantus[i][1]*SIZE)
lelantus.append(['Values','',8*T])
lelantus.append(['Fee','',8])

print('Data for Lelantus mint transaction')
print('T = {}'.format(T))
print()
print(tabulate.tabulate(lelantus, headers=headers))
total = sum([row[2] for row in lelantus])
print()
print('Total size is {} kB'.format(round(total/1024,2)))

print()
print()

# Triptych
#
# A mint transaction has the following components, with number of group/scalar elements
# (*) Note: fee and values are represented by 8-byte integers
#
# Output keys:				T
# Diffie-Hellman key:		1
# Values (*):				{8 bytes}*T
# Fee (*):					{8 bytes}
triptych = []
triptych.append(['Output keys', T])
triptych.append(['Diffie-Hellman key', 1])

for i in range(len(triptych)):
	triptych[i].append(triptych[i][1]*SIZE)
triptych.append(['Values','',8*T])
triptych.append(['Fee','',8])

print('Data for Triptych mint transaction')
print('T = {}'.format(T))
print()
print(tabulate.tabulate(triptych, headers=headers))
total = sum([row[2] for row in triptych])
print()
print('Total size is {} kB'.format(round(total/1024,2)))