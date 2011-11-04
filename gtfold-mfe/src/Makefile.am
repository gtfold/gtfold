INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

AM_CFLAGS = $(OPENMP_CFLAGS) -DDATADIR='$(datadir)/@PACKAGE@'
AM_CXXFLAGS = -DDATADIR='$(datadir)/@PACKAGE@'

bin_PROGRAMS = gtfold

gtfold_SOURCES = \
	main.cc\
	loader.cc\
	utils.cc\
	options.cc\
	constraints.cc\
	global.c\
	energy.c\
	algorithms.c\
	traceback.c\
	subopt_traceback.cc\
	stochastic-sampling.cc\
  algorithms-partition.c\
	partition-dangle.c\
	partition-func.c\
	shapereader.cc

gtfold_LDFLAGS = 

gtfold_LDADD = -lm

CLEANFILES = *~ *.o
