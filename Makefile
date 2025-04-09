CC = g++

ROOT_FLAGS = `root-config --cflags`
ROOT_LIBS = `root-config --libs`

SMDR_LIBS = -lm -lsmdr -ltsil -l3vil

all: init_random generate_random generate_net fit theor_err

init_random: init_random.cpp loop_Fit_Inputs.cpp loop_configs.cpp smdr_pdg_2025.h
	$(CC) $(ROOT_FLAGS) init_random.cpp $(ROOT_LIBS) $(SMDR_LIBS) -o init_random
	
generate_random: generate_random.cpp loop_Fit_Inputs.cpp loop_configs.cpp smdr_pdg_2025.h
	$(CC) $(ROOT_FLAGS) generate_random.cpp $(ROOT_LIBS) $(SMDR_LIBS) -o generate_random
	
generate_net: generate_net.cpp loop_Fit_Inputs.cpp loop_configs.cpp smdr_pdg_2025.h
	$(CC) $(ROOT_FLAGS) generate_net.cpp $(ROOT_LIBS) $(SMDR_LIBS) -o generate_net
	
fit: fit.cpp smdr_pdg_2025.h
	$(CC) $(ROOT_FLAGS) fit.cpp $(ROOT_LIBS) -o fit
	
theor_err: theor_err.cpp loop_Fit_Inputs.cpp loop_configs.cpp smdr_pdg_2025.h
	$(CC) theor_err.cpp $(SMDR_LIBS) -o theor_err
	
clean: 
	rm -rf init_random generate_random generate_net fit theor_err
