
help :
	@echo "make <target name>"
test1: test1.cpp ../utils/PRgenerator.cpp 
	g++ -g $^ -o target/$@ -Wno-pointer-arith -fopenmp

test2: test2.cpp ../utils/SocketMessenger.cpp
	g++ -g $^ -o target/$@ -Wno-pointer-arith -fopenmp

test_securenn: test_securenn.cpp ../utils/SocketMessenger.cpp ../utils/PRgenerator.cpp ../mpc_engine/SecurennEngine.cpp
	g++ -g $^ -o target/$@ -Wno-pointer-arith -fopenmp

test_ABY3: test_ABY3.cpp ../utils/SocketMessenger.cpp ../utils/PRgenerator.cpp ../mpc_engine/ABY3Engine.cpp
	g++ -g $^ -o target/$@ -Wno-pointer-arith -fopenmp

test_Myengine: test_Myengine.cpp ../utils/SocketMessenger.cpp ../utils/PRgenerator.cpp ../mpc_engine/MyEngine.cpp
	g++ -g $^ -o target/$@ -Wno-pointer-arith -fopenmp

run-test2: test2
	@./target/test2 0 > ./log/party0.txt &
	@./target/test2 1 > ./log/party1.txt &
	@./target/test2 2 > ./log/party2.txt 
	@echo "complete run-test2"

run-test_securenn: test_securenn
	@./target/test_securenn 0 > ./log/party0.txt &
	@./target/test_securenn 1 > ./log/party1.txt &
	@./target/test_securenn 2 > ./log/party2.txt 
	@echo "complete run-test_securenn"

run-test_ABY3: test_ABY3
	@./target/test_ABY3 0 > ./log/party0.txt &
	@./target/test_ABY3 1 > ./log/party1.txt &
	@./target/test_ABY3 2 > ./log/party2.txt 
	@echo "complete run-test_ABY3"

run-test_Myengine: test_Myengine
	@./target/test_Myengine 0 > ./log/party0.txt &
	@./target/test_Myengine 1 > ./log/party1.txt &
	@./target/test_Myengine 2 > ./log/party2.txt 
	@echo "complete run-test_Myengine"