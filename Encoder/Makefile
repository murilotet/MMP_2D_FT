CXX = g++
CXXFLAGS = -g -fopenmp -Wall -O3 -I ../CODARITMETICO/LONGLONG/ -I ../DICIONARIO/
CPPFLAGS = 

LDFLAGS = 

LDLIBS = -lm 

SOURCES = ../CODARITMETICO/LONGLONG/Aritmetico.cpp ../DICIONARIO/Dicionario2D-FT.cpp ./MMP2D-FT-MP.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: compile clean all

compile: MMP2D-FT-MP 
clean:
	$(RM) $(RMFLAGS) $(OBJECTS) MMP2D-FT-MP
all: clean compile
 
MMP2D-FT-MP: LINK.o = $(LINK.cc)
MMP2D-FT-MP: $(OBJECTS)
