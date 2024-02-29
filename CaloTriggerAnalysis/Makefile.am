AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

pkginclude_HEADERS = \
  MBDEmulatorTreeMaker.h \
  MBDEmulatorTreeMakerLinkDef.h \
  HCALEmulatorTreeMaker.h \
  HCALEmulatorTreeMakerLinkDef.h \
  CaloEmulatorTreeMakerLinkDef.h \
  CaloEmulatorTreeMaker.h \
  TrigTreeMaker.h \
  TrigTreeMakerLinkDef.h \
  WaveTreeMaker.h \
  WaveTreeMakerLinkDef.h


ROOTDICTS = 

pcmdir = $(libdir)
nobase_dist_pcm_DATA = 

lib_LTLIBRARIES = \
  libtrigtreemaker_io.la \
  libtrigtreemaker.la \
  libwavetreemaker.la \
  libmbdemulatortreemaker.la \
  libhcalemulatortreemaker.la \
  libcaloemulatortreemaker.la

libtrigtreemaker_io_la_SOURCES = \
  $(ROOTDICTS) 

libtrigtreemaker_la_SOURCES = \
  TrigTreeMaker.cc

libtrigtreemaker_io_la_LIBADD = \
  -lphool \
  -lll1out

libtrigtreemaker_la_LIBADD = \
  libtrigtreemaker_io.la \
  -lcalowaveformsim \
  -lll1out \
  -lcalo_io \
  -lphhepmc_io \
  -lffaobjects \
  -ljetbackground_io \
  -lgsl \
  -lgslcblas \
  -lphg4hit \
  -lcentrality_io \
  -lSubsysReco

libwavetreemaker_la_SOURCES = \
  WaveTreeMaker.cc

libwavetreemaker_la_LIBADD = \
  -lcalowaveformsim \
  -lcalo_io \
  -lffaobjects \
  -lSubsysReco


libmbdemulatortreemaker_la_SOURCES = \
  MBDEmulatorTreeMaker.cc

libmbdemulatortreemaker_la_LIBADD = \
  -lcalowaveformsim \
  -lll1out \
  -lSubsysReco

libhcalemulatortreemaker_la_SOURCES = \
  HCALEmulatorTreeMaker.cc

libhcalemulatortreemaker_la_LIBADD = \
  -lcalowaveformsim \
  -lll1out \
  -lSubsysReco

libcaloemulatortreemaker_la_SOURCES = \
  CaloEmulatorTreeMaker.cc

libcaloemulatortreemaker_la_LIBADD = \
  -lcalowaveformsim \
  -lll1out \
  -lSubsysReco

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals_io \
  testexternals \
  testexternals_wave \
  testexternals_mbd \
  testexternals_hcal \
  testexternals_calo

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD   = libtrigtreemaker_io.la

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libtrigtreemaker.la

testexternals_wave_SOURCES = testexternals.cc
testexternals_wave_LDADD   = libwavetreemaker.la

testexternals_mbd_SOURCES = testexternals.cc
testexternals_mbd_LDADD   = libmbdemulatortreemaker.la

testexternals_hcal_SOURCES = testexternals.cc
testexternals_hcal_LDADD   = libhcalemulatortreemaker.la

testexternals_calo_SOURCES = testexternals.cc
testexternals_calo_LDADD   = libcaloemulatortreemaker.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
