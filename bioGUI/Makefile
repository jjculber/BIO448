ALL = BioGUI

CC=javac
ENGINE=java
JFLAGS=-deprecation #-Xlint

all: $(ALL)

BioGUI: BioGUI.java
	@echo "-------------------------------"
	$(CC) $(JFLAGS) bioGUI/$^
	@echo "-------------------------------"

run:
	@echo "running..."
	$(ENGINE) bioGUI.BioGUI
	@echo "finished running."

BioGUI.java:

clean:
	@echo "-------------------------------"
	@echo "*** Cleaning Files..."
	rm -f bioGUI/*.class
	rm -f bioGUI/*/*.class
	@echo "-------------------------------"
