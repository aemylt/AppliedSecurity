target       : $(wildcard target.[ch])
	@gcc -std=gnu99 -o ${@} $(filter %.c, ${^})

attack       : $(wildcard attack.[ch])
	@gcc -std=gnu99 -o ${@} $(filter %.c, ${^})

attack.class : attack.java
	@javac ${^}

all   :        target attack attack.class

clean :
	@rm -f target attack attack.class *.pyc *.pyo
