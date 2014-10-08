
all: hapsplit

hapsplit: hapsplit.cpp GFA.h GFA.cpp
	g++ -o $@ $^
