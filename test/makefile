
help :
	@echo "make <target name>"
test1: test1.cpp ../utils/PRgenerator.cpp 
	g++ -g $^ -o target/$@ -Wno-pointer-arith -fopenmp

test2: test2.cpp ../utils/SocketMessenger.cpp
	g++ -g $^ -o target/$@ -Wno-pointer-arith -fopenmp

run-test2: test2
	@./target/test2 0 > ./log/party0.txt &
	@./target/test2 1 > ./log/party1.txt &
	@./target/test2 2 > ./log/party2.txt 
	@echo "complete run-test2"