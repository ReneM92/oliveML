# --- MACROS

# define name of executable
EXEC = olive 

# define the C/C++ compiler to use
CC= g++

# define any compile-time flags
CFLAGS= #-O3
LDFLAGS =


# define C/C++ source files
SRCS= sim.cpp axon.cpp caModel.cpp cell.cpp dendrite.cpp helper.cpp soma.cpp

# define C/C++ header files
HDRS= axon.hpp caModel.hpp cell.hpp dendrite.hpp helper.hpp soma.hpp

# --- TARGETS
all: ${EXEC}

#Builds the program
${EXEC}: ${SRCS} ${HDRS}
	@echo #
	@echo "-- BUILDING PROGRAM --"
	${CC} ${SRCS} ${CFLAGS} ${LDFLAGS}  -o ${EXEC}

clean:
	@echo #
	@echo "-- CLEANING PROJECT FILES --"
	$(RM) *.o ${EXEC} *.dat *.png

run: 
	./${EXEC}
