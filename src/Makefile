# http://www.rapidtables.com/code/linux/gcc/gcc-g.htm
# -g    : No debug info.
# -lm   : Link with math library.
# -ansi : the same as -std=c89



# Variable definitions.
CC = gcc
# CFLAGS = -Wall -std=c99 -pedantic -g
CFLAGS = -Wall -std=c99 -O3
# CFLAGS = -Wall -std=c99 -pedantic -g
BUILDDIR = ../bin
RESULTSDIR = ../results


.PHONY: all
all: mkdirs clean compile
	@echo MAKE Complete...


.PHONY: compile
compile: heatbugs.c heatbugs.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir $(BUILDDIR); fi
	$(CC) heatbugs.c heatbugs.h $(CFLAGS) `pkg-config --cflags --libs glib-2.0` -o $(BUILDDIR)/heatbugs


.PHONY: mkdirs
mkdirs:
#	@if [ ! -d $(BUILDDIR) ]; then mkdir -p $(BUILDDIR); fi
#	@if [ ! -d $(RESULTSDIR) ]; then mkdir -p $(RESULTSDIR); fi
	mkdir -p $(BUILDDIR)
	mkdir -p $(RESULTSDIR)


.PHONY: clean
clean:
	rm -rf $(BUILDDIR)/*
#	rm -d $(BUILDDIR)
#	rm -drf $(BUILDDIR)
