
all: gfademo

gfademo: gfademo.cpp GFA.h GFA.cpp
	g++ -o $@ $^
