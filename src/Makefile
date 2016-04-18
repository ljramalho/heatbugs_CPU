# Define required directories.
#export OBJDIR := ${CURDIR}/obj
#export BUILDDIR := ${CURDIR}/build
#SRCDIRS = ${CURDIR}/src

# -lm : Link with math library-


# Variable definitions.
CC = gcc
CFLAGS = -Wall -std=c99 -pedantic -g
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
