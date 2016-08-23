# -lm : Link with math library.

# http://www.rapidtables.com/code/linux/gcc/gcc-g.htm
# -g :  No debug info.

# -ansi : the same as -std=c89



# Variable definitions.
CC = gcc
# CFLAGS = -std=c99 -Wall -pedantic -g
CFLASG = -std=C99 -Wall -g
# CFLASG = -std=C99 -Wall -O3
BUILDDIR = ../bin


.PHONY: all
all: clean compile
	@echo MAKE Complete...

.PHONY: compile
compile: heatbugs.c heatbugs.h
	@if [ ! -d $(BUILDDIR) ]; then mkdir $(BUILDDIR); fi
	$(CC) heatbugs.c heatbugs.h $(CFLAGS) `pkg-config --cflags --libs glib-2.0` -o $(BUILDDIR)/heatbugs
#	$(CC) heatbugs.c heatbugs.h $(CFLAGS) `pkg-config --cflags --libs cf4ocl2 glib-2.0` -lOpenCL -o $(BUILDDIR)/heatbugs


.PHONY: mkdirs
mkdirs:
	mkdir -p $(OBJDIR)
	mkdir -p $(BUILDDIR)

.PHONY: clean
clean:
	rm -rf $(BUILDDIR)/heatbugs
#	rm -d $(BUILDDIR)
#	rm -drf $(BUILDDIR)
