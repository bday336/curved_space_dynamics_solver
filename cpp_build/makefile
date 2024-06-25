
# Define the symbols we might want to change:
CXX	:= g++
CXXFLAGS := -std=c++20 -lginac -lcln
INCLUDES := /includes/eigen-3.4.0/Eigen/		# This would contain any -I options to the compiler, if there are any.

OBJECTS	:= $(patsubst %.cpp,%.o,$(wildcard *.cpp Computation/*.cpp Computation/Solvers/*.cpp AmbientSpace/*.cpp AmbientSpace/Components/*.cpp AmbientSpace/ExampleSpaces/*.cpp ConfigurationSpace/*.cpp))

simulation_test: $(OBJECTS)
	$(CXX) $(CXXFLAGS) -I $(INCLUDES)  $(OBJECTS) -o simulation_test 

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -I $(INCLUDES) -c $< -o $@

clean:
	rm $(OBJECTS) simulation_test	









# objects/%.o: Computation/%.cpp
# 	$(CXX) $(CXXFLAGS) -c $<

# State.o: State.cpp
# 	g++ -c State.cpp -std=c++11

# DState.o: DState.cpp
# 	g++ -c DState.cpp -std=c++11

# DataList.o: DataList.cpp
# 	g++ -c DataList.cpp -std=c++11

# my_program: $(OBJECTS)
# 	$(CXX) $(OBJECTS) -o my_program

# processing.o: processing.cxx
# 	$(CXX) $(INCLUDES) $(CXXFLAGS) -c processing.cxx -o processing.o

# gui.o: gui.cxx
# 	$(CXX) $(CXXFLAGS) -c gui.cxx -o gui.o




# simulation_test: run.o Computation/State.o Computation/DState.o Computation/DataList.o
# 	g++ run.o Computation/State.o Computation/DState.o Computation/DataList.o -o simulation_test -std=c++11

# run.o: run.cpp
# 	g++ -c run.cpp -std=c++11

# State.o: State.cpp
# 	g++ -c State.cpp -std=c++11

# DState.o: DState.cpp
# 	g++ -c DState.cpp -std=c++11

# DataList.o: DataList.cpp
# 	g++ -c DataList.cpp -std=c++11


# clean:
# 	rm *.o Computation/*.o simulation_test	