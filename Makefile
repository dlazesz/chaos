INSTALL_DIR = $(HOME)/vnodelp_bin
VNODE_DIR = $(HOME)/vnodelp
CONFIG_FILE = $(VNODE_DIR)/config/LinuxWithProfil

include $(CONFIG_FILE)

CXXFLAGS += -I$(INSTALL_DIR)/include\
	-I$(INSTALL_DIR)/FADBAD++
LDFLAGS  += -L$(INSTALL_DIR)/lib

all: chaos_GSL chaos_P00 chaos_P83

chaos_GSL:
	g++ chaos_GSL.cc -lgsl -lgslcblas -o enso

chaos_P00:
	chaos_P00: chaos_P00.o
	$(CXX) $(LDFLAGS) -o $@ chaos_P00.o -lvnode $(LIBS)

chaos_P83:	
	chaos_P83: chaos_P83.o
	$(CXX) $(LDFLAGS) -o $@ chaos_P83.o -lvnode $(LIBS)

clean:
	rm -rf chaos_GSL chaos_P00 chaos_P83 *.o
