bin_PROGRAMS = varval

varval_CPPFLAGS = \
     -I$(top_srcdir)/SeqLib/src \
     -I$(top_srcdir)/SeqLib/htslib -Wno-sign-compare \
     -I$(top_srcdir)/SeqLib 

varval_LDADD = \
	$(top_builddir)/SeqLib/src/libseqlib.a \
	$(top_builddir)/SeqLib/bwa/libbwa.a \
	$(top_builddir)/SeqLib/htslib/libhts.a -llzma -lbz2

##varval_SOURCES = DiscordantCluster.cpp AlignedContig.cpp BreakPoint.cpp \
##			varval.cpp svabaRead.cpp STCoverage.cpp AlignmentFragement.cpp \
##			gzstream.C

varval_SOURCES = varval.cpp gzstream.C ValidatorBamReader.cpp
