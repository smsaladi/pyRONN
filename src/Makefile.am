bin_PROGRAMS = RONN

#RONN application
RONN_SOURCES = RONN.cpp callBBF.cpp
RONN_CPPFLAGS = -I./include -O2 -fopenmp
RONN_LDADD = -lm -lgomp

install-exec-hook:
	ln -s $(prefix)/bin/RONN /usr/local/bin/RONN_32
	cp -R ../data/* $(prefix)

uninstall-hook:
	rm -rf $(prefix)
	rm /usr/local/bin/RONN_32
