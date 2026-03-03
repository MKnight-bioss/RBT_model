#ifndef SIM_MANAGER_H
#define SIM_MANAGER_H

#include <string>
#include <unordered_map>

struct SimResult;

struct SimManager{
	
	SimManager(const std::string& paramFile, const std::string& farmDataFile);
	
	SimResult runSingleSim(int sim_id, const std::unordered_map<std::string, double>& overrides);
	void runSims(int numSims, const std::string& outputFile, const std::unordered_map<std::string, double>& overrides = {});
	
	std::string param_file;
	std::string farm_data_file;
};

#endif