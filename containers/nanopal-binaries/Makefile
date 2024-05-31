.PHONY = all

all: bin/cigar_parser bin/inter bin/RM_collapse

bin/cigar_parser: cigar_parser.cpp Makefile
	mkdir -p bin
	g++ -o $@ -Wall -O2 $<

bin/inter: inter.0118.cpp Makefile
	mkdir -p bin
	g++ -o $@ -Wall -O2 $<

bin/RM_collapse: RM_collapse.0921.cpp Makefile
	mkdir -p bin
	g++ -o $@ -Wall -O2 $<
