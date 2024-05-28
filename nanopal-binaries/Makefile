.PHONY = all

all: cigar_parser

cigar_parser: cigar_parser.cpp Makefile
	mkdir -p bin
	g++ -o bin/$@ -Wall -O2 $<
