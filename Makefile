#OBJS specifies which files to compile as part of the project
#OBJS = barnes_hut.cpp
OBJS = nbody-simulator.cpp barnes_hut.cpp

#CC specifies which compiler we're using
#CC = g++
CC = x86_64-w64-mingw32-g++

#Compiler flags
COMPILER_FLAGS = -Wall -pedantic

#LINKER_FLAGS specifies the libraries we're linking against
LINKER_FLAGS = -lmingw32 -lglu32 -lopengl32 -lglut

#OBJ_NAME specifies the name of our exectuable
#OBJ_NAME = barnes.exe
OBJ_NAME = nbody.exe

#This is the target that compiles our executable
all : $(OBJS) 
		$(CC) $(OBJS) $(LINKER_FLAGS) $(COMPILER_FLAGS) -o $(OBJ_NAME)