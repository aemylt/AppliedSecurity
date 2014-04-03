encrypt       : $(wildcard encrypt.[ch])
	@gcc -std=gnu99 -o ${@} $(filter %.c, ${^}) -lcrypto

encrypt.class : encrypt.java
	@javac ${^}

all   :        encrypt encrypt.class

clean :
	@rm -f encrypt encrypt.class *.pyc *.pyo
