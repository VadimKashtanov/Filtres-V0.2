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

#gdb 
#time
./prog

#valgrind --leak-check=yes --track-origins=yes ./prog