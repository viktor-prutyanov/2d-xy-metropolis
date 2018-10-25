CXXFLAGS:=-Wall -O3 -std=c++14 -mkl -xHost
LDFLAGS:=-lm

src:=main.cpp

CXX:=icpc

all: xy

%.d: %.cpp
	$(CXX) -MM $(CXXFLAGS) -o $@ $<

-include $(src:.cpp=.d)

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $(<:.d=.cc)

xy: $(src:.cpp=.o)
	$(CXX) $^ $(CXXFLAGS) -o $@ $(LDFLAGS)

tags:
	rm -f tags
	find . -name '*.cpp' -exec ctags --append {} +

clean:
	rm -f *.o *.d xy

.PHONY: all clean tags
