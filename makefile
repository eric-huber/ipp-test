CC        = g++
CXXFLAGS += -I /opt/intel/compilers_and_libraries_2016.0.109/linux/ipp/include/
CXXFLAGS += -std=c++11
LDFLAGS  += -lboost_program_options
LDFLAGS  += -lm
LDFLAGS  += -L/opt/intel/compilers_and_libraries_2016.0.109/linux/ipp/lib/intel64 -lippcc -lipps -lippcore

PROG=ipp-test
OBJS=main.o

.PHONY: all clean
$(PROG): $(OBJS)
	$(CC) -o $(PROG) $(OBJS) $(LDFLAGS)

%.o: %.cc
	$(CC) -c $(CXXFLAGS) $<

all: $(PROG)

clean:
	rm -f $(OBJS) $(PROG) fft-*.txt
