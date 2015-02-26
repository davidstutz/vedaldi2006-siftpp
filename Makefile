# file:        Makefile
# author:      Andrea Vedaldi
# description: Build SIFT++

CXXFLAGS += -g -O3
DIST=siftpp
VER=0.2.2

# --------------------------------------------------------------------
#
# --------------------------------------------------------------------

# Determine on the flight the system we are running on
Darwin_ARCH := mac
Linux_ARCH  := glx
ARCH := $($(shell uname)_ARCH)

BINDIST=$(DIST)-$(ARCH)-$(VER)

# --------------------------------------------------------------------
#
# --------------------------------------------------------------------

sift : sift.o sift-driver.o
	g++ $^ -o $@

.PHONY: clean
clean:
	rm -f *.o
	find . -name '*~' -exec rm -f \{\} \;
	find . -name '.DS_Store' -exec rm -f \{\} \;

.PHONY: distclean
distclean: clean
	rm -rf results
	rm -f sift
	rm -f .gdb_history
	rm -f $(DIST)-*.tar.gz
	rm -rf $(BINDIST)
	rm -f $(BINDIST).tar.gz

.PHONY: test
test: sift
	(test -e results || mkdir results)
	./sift --verbose --output results/img3 --save-gss data/img3.pgm

.PHONY: autorights
autorights:
	autorights . \
	  --verbose \
	  --recursive \
	  --template notice.txt \
	  --years 2006 \
	  --authors "Andrea Vedaldi (UCLA VisionLab)" \
	  --program "SIFT++"


.PHONY: dist
dist: distclean
	echo Version $(VER) - Archived on `date` > TIMESTAMP
	tar chvzf $(DIST)-$(VER).tar.gz ../$(DIST)

.PHONY: bindist
bindist: sift
	test -e $(BINDIST) || mkdir $(BINDIST)
	cp sift $(BINDIST)
	cd $(BINDIST) ; strip -S sift
	tar cvzf $(BINDIST).tar.gz $(BINDIST)
