CPPFLAGS  := -g -O3 --std=c++17

FAKEROOT  := $(HOME)/fakeroot
INCLUDES  := -I$(FAKEROOT)/include
LIBRARIES := -L$(FAKEROOT)/lib

LIBS_HS107 := -lcasadi -lcasadi_nlpsol_ipopt
LIBS_GRI30 := -lcantera -lblas -lcblas -llapack -lpthread

all: hs107.cpp
	g++ $(CPPFLAGS) $(INCLUDES) hs107.cpp $(LIBRARIES) $(LIBS_HS107) -o hs107.exe
	g++ $(CPPFLAGS) $(INCLUDES) gri30.cpp $(LIBRARIES) $(LIBS_GRI30) -o gri30.exe

clean:
	rm -rf *.exe
