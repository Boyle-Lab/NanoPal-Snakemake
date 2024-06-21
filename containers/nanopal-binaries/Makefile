.PHONY = all

all: bin/cigar_parser bin/inter bin/RM_collapse bin/cluster bin/collapse bin/intersect

bin/cigar_parser: cigar_parser.cpp Makefile
	mkdir -p bin
	g++ -o $@ -Wall -O2 $<

bin/inter: inter.0118.cpp Makefile
	mkdir -p bin
	g++ -o $@ -Wall -O2 $<

bin/RM_collapse: RM_collapse.0921.cpp Makefile
	mkdir -p bin
	g++ -o $@ -Wall -O2 $<

bin/cluster: cluster.1219.cpp Makefile
	mkdir -p bin
	g++ -o $@ -Wall -O2 $<

bin/intersect: intersect.lisp Makefile
	mkdir -p bin
	sbcl --disable-debugger --load $< --eval '(intersect::build "'$@'")'

bin/collapse: collapse.lisp Makefile
	mkdir -p bin
	sbcl --disable-debugger --load $< --eval '(collapse::build "'$@'")'
