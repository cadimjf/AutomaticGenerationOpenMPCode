INC=	-I /usr/include/libxml2/

LIBS=	-lxml2

COMPILER= g++


all:	main.o
	$(COMPILER) $(INC) $(LIBS) -o agoscc main.o -w
	rm *.o

%.o:	%.cpp
	$(COMPILER) -O3  -c $(INC) -w $<
