
SRC=gfademo.cpp GFA.cpp
HEADER=GFA.h

all: gfademo

gfademo: $(SRC) $(HEADER)
	g++ -o $@ $(SRC)
