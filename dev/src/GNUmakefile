PROGNAME = ../bin/minimal_test
SOURCEFILES = main.cc sWeights.cc D2hhmumuReader.cc D2KKmumuReader.cc D2KpimumuReader.cc D2pipimumuReader.cc D2hhmumuFitter.cc RooJohnsonSU.cc RooThreshold.cc D2hhmumuModel.cc
OBJS    = $(patsubst %.cc, %.o,$(SOURCEFILES))

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

LDFLAGS       = -O
LIBS         += $(ROOTLIBS)
CFLAGS       += $(ROOTCFLAGS) -std=c++0x

#  Not sure why Minuit isn't being included -- put in by hand
#  #
LIBS         += -lMinuit -L /opt/root/lib \
			        -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix \
	   		        -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic \
				        -lSpectrum -lRooFit -lHtml -l RooStats 

%.o: %.cc
	        g++ ${CFLAGS} -c  -g -o $@ $<

$(PROGNAME):$(OBJS)
	        g++ -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
		#-mv  $(OBJS) ../obj

test:
	        @echo $(ROOTCFLAGS)

clean:  
	        -rm -f ${PROGNAME} ${OBJS}
