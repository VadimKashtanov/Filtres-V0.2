clear
printf "[***] Compilation ... \n"

#	Compiler
#gcc -g src/*.c src/selection/gen/*.c src/selection/sel/*.c -Ihead -o prog -lm -lSDL2 -lSDL2main -lSDL2_ttf #-fopenmp
python3 compiler.py
if [ $? -eq 1 ]
then
	printf "[***] Erreure. Pas d'execution.\n"
	exit
fi

#	Executer
printf "[\033[95m***\033[0m] ========= Execution du programme =========\n"

valgrind --leak-check=yes --track-origins=yes ./prog
#time ./prog
if [ $? -ne 0 ]
then
	printf "[***] Erreur durant l'execution.\n"
	#valgrind --leak-check=yes --track-origins=yes ./prog
	systemd-run --scope -p MemoryMax=100M gdb ./prog
	exit
fi