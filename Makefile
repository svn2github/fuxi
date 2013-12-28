SMITHLAB_CPP = ./smithlab_cpp

OPT = 1

ifndef SMITHLAB_CPP
$(error Must define SMITHLAB_CPP variable)
endif

PROGS = extract_multiple_run_kcv generate_hash_function build_metagenome_hash_table build_graph  angle_measurement \
	    build_graph_on_cluster visualize_metagenome exec_query db_insert angle_measurement simulate_fastq generate_metagenome_study \
	    mix_metagene

SOURCES = $(wildcard *.cpp)

INCLUDEDIRS = $(SMITHLAB_CPP)
INCLUDEARGS = $(addprefix -I,$(INCLUDEDIRS))

LIBS += -lgsl -lgslcblas -lboost_filesystem -lboost_system

CXX = g++
CXXFLAGS = -Wall -fmessage-length=50
OPTFLAGS = -O2
DEBUGFLAGS = -g

ifeq "$(shell uname)" "Darwin"
CFLAGS += -arch x86_64
endif

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

all: $(PROGS)

extract_kmer_counts : Metagenome.o

generate_hash_function : MGHashFunction.o MGFeatureSet.o Metagenome.o

build_metagenome_hash_table : MGHashTable.o MGHashFunction.o MGFeatureSet.o Metagenome.o

simulate_fastq : NGS_simulator.o

build_graph_on_cluster : MGraph.o Metagenome.o MGHashTable.o

build_graph : MGraph.o Metagenome.o MGHashTable.o

visualize_metagenome : Metagenome.o

exec_query : MGraph.o Metagenome.o MGHashTable.o MGHashFunction.o

db_insert : MGraph.o Metagenome.o MGHashTable.o MGHashFunction.o

angle_measurement : Metagenome.o

convert_metagenome : Metagenome.o

angle_measurement : Metagenome.o

mix_metagene : NGS_simulator.o Metagenome.o

$(PROGS): $(addprefix $(SMITHLAB_CPP)/, GenomicRegion.o \
	smithlab_os.o smithlab_utils.o OptionParser.o)

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDEARGS)

%: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ $^ $(INCLUDEARGS) $(LIBS)

clean:
	@-rm -f $(PROGS) *.o *.so *.a *~

.PHONY: clean
