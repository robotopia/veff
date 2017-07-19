INSTALL_DIR = /usr/local/bin

veff: veff.c
	gcc -Wall -Wextra -o $@ $< -lpsrcat -lm

install:
	cp veff $(INSTALL_DIR)

clean:
	$(RM) veff
