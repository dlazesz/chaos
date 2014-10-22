INSTALL_DIR = $(HOME)/vnodelp_bin
CONFIG_FILE = $(INSTALL_DIR)/vnodelp/config/LinuxWithProfil

include $(CONFIG_FILE)

CXXFLAGS += -I$(INSTALL_DIR)/vnodelp/include\
	-I$(INSTALL_DIR)/vnodelp/FADBAD++
LDFLAGS  += -L$(INSTALL_DIR)/vnodelp/lib 

all: chaos_GSL chaos_P00 chaos_P83 chaos_P00_Z0.3995

chaos_GSL:
	$(CXX) -std=c++11 chaos_GSL.cc -lstdc++ -lm -lgsl -lgslcblas -o $@

chaos_P00:	chaos_P00.o
	$(CXX) $(LDFLAGS) -o $@ chaos_P00.o -lvnode $(LIBS)

chaos_P83:	chaos_P83.o
	$(CXX) $(LDFLAGS) -o $@ chaos_P83.o -lvnode $(LIBS)

chaos_P00_Z0.3995:	chaos_P00_Z0.3995.o
	$(CXX) $(LDFLAGS) -o $@ chaos_P00_Z0.3995.o -lvnode $(LIBS)

clean:
	rm -rf chaos_GSL chaos_P00 chaos_P83 chaos_P00_Z0.3995 *.o
