from os import listdir
import json
from os import system

def lire_somme(chemin):
	with open(chemin, 'rb') as co:
		elm = co.read()
		return [(i+1)*elm[i] for i in range(len(elm))]

chemins = [
	'impl/', #'impl/couches/'
	#'impl/mdl/', #'impl/mdl/coeur/', 'impl/mdl/gradient/', 'impl/mdl/outils/',
	#'impl/selection/',
	#'impl/insts/',
	#'impl/optimisation/',
	#'impl/prospection/'
]	#, 'src/selection/gen/', 'src/selection/sel/'] #ou utiliser walk()
fichiers = [
	chemin+fichier
		for chemin in chemins
			for fichier in listdir(chemin)
				if fichier[-2:]=='.c'
]

sommes = {fichier : lire_somme(fichier) for fichier in fichiers}
sommes_tete = {'def/'+tete:lire_somme('def/'+tete) for tete in listdir('def')}

a_compiler = []

if not 'trace_tete' in listdir():
	#	Pas de trace de la tete, on l'ecrit et on compile tout
	a_compiler = fichiers
	with open('trace_tete', 'w') as co: co.write(json.dumps(sommes_tete))
	with open('trace', 'w') as co: co.write(json.dumps(sommes))
else:
	with open('trace_tete', 'r') as co: sommes_de_trace_tete = eval(co.read())
	#print(sommes_tete, sommes_de_trace_tete)
	
	if not all('def/'+tete in sommes_de_trace_tete.keys() for tete in listdir('def')):
		a_compiler = fichiers
	elif not all(sommes_de_trace_tete['def/'+tete] == sommes_tete['def/'+tete] for tete in listdir('def')):
		a_compiler = fichiers
	elif not 'trace' in listdir():
		#	Pas de fichier de trace on compile tout et on ecrit la trace
		with open('trace', 'w') as co: co.write(json.dumps(sommes))
		a_compiler = fichiers
	else:
		# Fichier trace -> on recompile que les fichier modifiees
		with open('trace', 'r') as co: sommes_de_trace = eval(co.read())
		for fichier,_somme in sommes.items():
			if not fichier in sommes_de_trace.keys():
				a_compiler += [fichier]
			elif sommes_de_trace[fichier] != _somme:
				a_compiler += [fichier]
			#else:
			#	on recompile pas

if not 'prog' in listdir():
	a_compiler = sommes

SDL_FLAGS = " -g " #"-O3" #" -mavx -g" #"-lSDL2 -lSDL2main -lSDL2_ttf `sdl2-config --cflags --libs`"

print(f"[***] {len(a_compiler)}/{len(sommes)} fichiers a compiler")
if len(a_compiler) == 0:
	exit(0)
for _a_compiler in a_compiler:
	nom_fichier = _a_compiler.split('/')[-1]
	e = system(f"gcc -c {_a_compiler} -Idef -lm" + SDL_FLAGS)
	if e != 0: exit(1)
	system(f"mv {nom_fichier.replace('.c','.o')} bin/{_a_compiler.replace('/','_').replace('.c','.o')}")

system("rm prog")
system("gcc bin/*.o -o prog -Idef -lm " + SDL_FLAGS)

##	Mise a jour si tout se passe bien
sommes = {fichier : lire_somme(fichier) for fichier in fichiers}
sommes_tete = {'def/'+tete:lire_somme('def/'+tete) for tete in listdir('def')}
with open('trace_tete', 'w') as co: co.write(json.dumps(sommes_tete))
with open('trace', 'w') as co: co.write(json.dumps(sommes))