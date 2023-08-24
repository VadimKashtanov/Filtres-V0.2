#! /usr/bin/python3

def normaliser(lst):
	_max = max(lst)
	_min = min(lst)
	return [(i-_min)/(_max-_min) for i in lst], _max, _min

def lire(fichier, col=2):
	with open(fichier, 'r') as co:
		lignes = co.read().replace('\n\n','').split('\n')[1:]
		return [float(ligne.split('\t')[col]) for ligne in lignes]

from sys import argv

#'EURUSD_H1.csv'
prixs, _max, _min = lire(argv[1]), None, None

import struct as st

#"prixs.bin"
with open(argv[2], "wb") as co:
	#co.write(st.pack( 'ffI', _max, _min, len(prixs) ))
	co.write(st.pack('I', len(prixs)))
	co.write(st.pack('f'*len(prixs), *prixs))
