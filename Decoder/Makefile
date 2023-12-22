CXX = g++
CXXFLAGS = -g -fopenmp -Wall -O3 -I ../CODARITMETICO/LONGLONG/ -I ../DICIONARIO/
CPPFLAGS = 

LDFLAGS = 

LDLIBS = -lm 

SOURCES = ../CODARITMETICO/LONGLONG/Aritmetico.cpp ../DICIONARIO/Dicionario2D-FT.cpp ./UNMMP2D-FT-MP.cpp
OBJECTS = $(SOURCES:.cpp=.o)

.PHONY: compile clean all

compile: UNMMP2D-FT-MP 
clean:
	$(RM) $(RMFLAGS) $(OBJECTS) UNMMP2D-FT-MP
all: clean compile
 
UNMMP2D-FT-MP: LINK.o = $(LINK.cc)
UNMMP2D-FT-MP: $(OBJECTS)
