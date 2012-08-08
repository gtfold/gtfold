INCLUDES = -I$(top_srcdir) -I$(top_srcdir)/include

AM_CFLAGS = $(OPENMP_CFLAGS) -DDATADIR='$(datadir)/@PACKAGE@'
AM_CXXFLAGS = -DDATADIR='$(datadir)/@PACKAGE@'

bin_PROGRAMS = gtfold

gtfold_SOURCES = \
	main.cc\
	mfe_main.cc\
  loader.cc\
	utils.cc\
	options.cc\
	constraints.cc\
	global.cc\
	energy.c\
	algorithms.c\
	traceback.c\
	subopt_main.cc\
	subopt_traceback.cc\
	pf-shel-check.cc\
	key.cc\
	stochastic-sampling.cc\
	stochastic-sampling-d2.cc\
  algorithms-partition.c\
  boltzmann_main.cc\
	partition-dangle.c\
	partition-func.c\
	partition-func-d2.cc\
	shapereader.cc\
gtfold_LDFLAGS = 

gtfold_LDADD = -lm

CLEANFILES = *~ *.o
