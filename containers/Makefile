.phony: all

all: minimap2.sif palmer.sif nanopal-binaries.sif samtools.sif liger2liger.sif minimera.sif

minimap2.sif: minimap2.def
	singularity build --fakeroot $@ $<

palmer.sif: palmer.def
	singularity build --fakeroot $@ $<

nanopal-binaries.sif: nanopal-binaries.def
	singularity build --fakeroot $@ $<

samtools.sif: samtools.def
	singularity build --fakeroot $@ $<

minimera.sif: minimera.def minimera/minimera minimera/minimera.1
	singularity build --fakeroot $@ $<

liger2liger.zip:
	wget https://github.com/rlorigro/Liger2LiGer/archive/refs/heads/main.zip -O liger2liger.zip

autoconf-2.69.tar.gz:
	curl -O http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz

liger2liger.tar.gz: liger2liger.zip
	rm -rf ./Liger2Liger-main/
	unzip liger2liger.zip
	tar -czf liger2liger.tar.gz Liger2LiGer-main/

liger2liger.sif: liger2liger.def minimap2.sif liger2liger.tar.gz autoconf-2.69.tar.gz liger2liger.sh
	singularity build --fakeroot $@ $<

blast.sif:
	singularity pull docker://staphb/blast:2.15.0

