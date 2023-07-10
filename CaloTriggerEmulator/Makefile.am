AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
	-I$(includedir) \
	-I$(OFFLINE_MAIN)/include \
	-I$(ROOTSYS)/include

lib_LTLIBRARIES = \
	libll1out.la \
	libmbdtriggeremulator.la \
	libmbdpmtbuilder.la \
	libcalotriggeremulator.la \
	libcalopacketgetter.la

AM_LDFLAGS = -L$(libdir) -L$(OFFLINE_MAIN)/lib

libll1out_la_LIBADD = \
	-lphool \
	-lfun4all

libmbdtriggeremulator_la_LIBADD = \
	-lphool \
	-lfun4all \
	-lg4detectors \
	-lonnxruntime \
	-lcalowaveformsim \
	libll1out.la

libmbdpmtbuilder_la_LIBADD = \
	-lphool \
	-lfun4all \
	-lg4detectors \
	-lonnxruntime \
	-lcalowaveformsim \
	libll1out.la

libcalotriggeremulator_la_LIBADD = \
	-lphool \
	-lfun4all \
	-lg4detectors \
	-lonnxruntime \
	libll1out.la \
	-lcalowaveformsim

libcalopacketgetter_la_LIBADD = \
	-lphool \
	-lfun4all \
	-lg4detectors \
	-lonnxruntime \
	-lcalowaveformsim

pkginclude_HEADERS = \
	LL1Defs.h \
	LL1Out.h \
	LL1Outv1.h \
	LL1Outv2.h \
	MBDTriggerEmulator.h \
	CaloTriggerEmulator.h \
	MBDPmtBuilder.h \
	CaloPacketGetter.h \
	TriggerDefs.h \
	TriggerPrimitive.h \
	TriggerPrimitiveContainer.h \
	TriggerPrimitiveContainerv1.h

ROOT_DICTS = \
	LL1Out_Dict.cc \
	LL1Outv1_Dict.cc \
	LL1Outv2_Dict.cc \
	TriggerPrimitive_Dict.cc \
	TriggerPrimitiveContainer_Dict.cc \
	TriggerPrimitiveContainerv1_Dict.cc

pcmdir = $(libdir)

nobase_dist_pcm_DATA = \
	LL1Out_Dict_rdict.pcm \
	LL1Outv1_Dict_rdict.pcm \
	LL1Outv2_Dict_rdict.pcm \
	TriggerPrimitive_Dict_rdict.pcm \
	TriggerPrimitiveContainer_Dict_rdict.pcm \
	TriggerPrimitiveContainerv1_Dict_rdict.pcm


libll1out_la_SOURCES = \
	$(ROOT_DICTS) \
	LL1Out.cc \
	LL1Outv1.cc \
	LL1Outv2.cc \
	TriggerDefs.cc \
	TriggerPrimitive.cc \
	TriggerPrimitiveContainer.cc \
	TriggerPrimitiveContainerv1.cc 

libcalotriggeremulator_la_SOURCES = \
	CaloTriggerEmulator.cc

libmbdtriggeremulator_la_SOURCES = \
	MBDTriggerEmulator.cc

libmbdpmtbuilder_la_SOURCES = \
	MBDPmtBuilder.cc

libcalopacketgetter_la_SOURCES = \
	CaloPacketGetter.cc

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

%_Dict_rdict.pcm: %_Dict.cc ;


################################################
# linking tests

noinst_PROGRAMS = \
	testexternals_calotriggeremulator \
	testexternals_calopacketgetter \
	testexternals_mbdtriggeremulator \
	testexternals_mbdpmtbuilder \
	testexternals_ll1out

BUILT_SOURCES = testexternals.C

testexternals_calotriggeremulator_SOURCES = testexternals.C
testexternals_calotriggeremulator_LDADD = libcalotriggeremulator.la

testexternals_calopacketgetter_SOURCES = testexternals.C
testexternals_calopacketgetter_LDADD = libcalopacketgetter.la

testexternals_mbdtriggeremulator_SOURCES = testexternals.C
testexternals_mbdtriggeremulator_LDADD = libmbdtriggeremulator.la

testexternals_mbdpmtbuilder_SOURCES = testexternals.C
testexternals_mbdpmtbuilder_LDADD = libmbdpmtbuilder.la

testexternals_ll1out_SOURCES = testexternals.C
testexternals_ll1out_LDADD = libll1out.la


testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
