#CC			  = /opt/intel/bin/icc
CC        = g++
CXXFLAGS += -O2 -std=c++11
LDFLAGS  += -lboost_program_options
LDFLAGS  += -lm
LDFLAGS  += -lippcc -lipps -lippcore

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
