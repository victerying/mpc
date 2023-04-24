UTILS_CPP = $(wildcard ../utils/*.cpp)
MPC_CPP = $(wildcard ../mpc_engine/*.cpp)

help:
	@echo "this is a help message."
	@echo "UTILS_CPP: $(UTILS_CPP)"
	@echo "MPC_CPP: $(MPC_CPP)"


benchmark1: benchmark1.cpp $(UTILS_CPP) $(MPC_CPP)
	g++ -g $^ -o target/$@ -Wno-pointer-arith -fopenmp

run-benchmark1:
	@./target/benchmark1 0 > ./log/party0.txt &
	@./target/benchmark1 1 > ./log/party1.txt &
	@./target/benchmark1 2 > ./log/party2.txt 
	@echo "complete run-benchmark1"