AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
	libcalowaveformsim.la

AM_LDFLAGS = \
	-L$(libdir) \
	-L$(OFFLINE_MAIN)/lib

AM_CPPFLAGS = \
	-I$(includedir) \
	-I$(OFFLINE_MAIN)/include \
	-I$(ROOTSYS)/include

libcalowaveformsim_la_LIBADD = \
	-L$(libdir) \
	-L$(OFFLINE_MAIN)/lib \
	-lcalo_io \
	-lphool \
	-lfun4all \
	-lg4detectors_io \
	-lg4detectors \
	-lsph_onnx \
	-lcalo_reco \
	-lmbd_io \
	-lphg4hit

pkginclude_HEADERS = \
	ADCDefs.h \
	WaveformReturnCodes.h \
	WaveformContainer.h \
	WaveformContainerv1.h \
	CaloWaveFormSim.h \
	CaloWaveFormToy.h

pkgdata_DATA = \
	no_raddmgnoise.csv \
	low_raddmgnoise.csv \
	medium_raddmgnoise.csv \ 
	high_raddmgnoise.csv \
	testbeam_emcal_template.root \
	testbeam_ihcal_template.root \
	testbeam_ohcal_template.root \
	laser_mbd_template.root 

pcmdir = $(libdir)

nobase_dist_pcm_DATA = \
	WaveformContainer_Dict_rdict.pcm \
	WaveformContainerv1_Dict_rdict.pcm

libcalowaveformsim_la_SOURCES = \
	WaveformContainer_Dict.cc \
	WaveformContainerv1_Dict.cc \
	WaveformContainer.cc \
	WaveformContainerv1.cc \
	CaloWaveFormSim.cc \
	CaloWaveFormToy.cc


%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
	testexternals_calowaveformsim

testexternals_calowaveformsim_SOURCES = testexternals.C
testexternals_calowaveformsim_LDADD = libcalowaveformsim.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating table CINT dictionaries.
clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
