#include "SimManager.h"
#include <iostream>
#include <unordered_map>
#include <vector>
#include <algorithm>

std::string makeSafeFilename(double val){
	
	std::string s = std::to_string(val);
	s.erase(s.find_last_not_of('0') + 1, std::string::npos);
	
	if (!s.empty() && s.back() == '.')
		s.pop_back();
	
	size_t pos = s.find('.');
	if (pos != std::string::npos)
		s.replace(pos,1,"pt");
	
	return s;
}

int main(){
	
	int num_params = 21;
	int num_reps = 10;
	
	std::vector<double> frac_not_sharing;
	frac_not_sharing.resize(num_params, 0.0);
	
	for (std::size_t i = 0; i < frac_not_sharing.size(); ++i){
		frac_not_sharing[i] = 0.05*i;
	}
	
	for (int i = 0; i < num_params; ++i){
		
		auto time_start = time(NULL);
		
		std::cout << "Starting parameter " << i+1 << std::endl;
		
		std::unordered_map<std::string, double> overrides;
		overrides["frac_sellers_not_sharing"] = frac_not_sharing[i];
		overrides["risk_score_not_sharing"] = 0.0;
		
		std::string param_file = "sim_params.csv";
		std::string farm_data_file = "farm_data.csv";
		std::string output_file = "avg_time_series_output_frac_not_sharing_" + makeSafeFileName(frac_not_sharing[i]) + "_risk_score_not_sharing_0.csv";
		
		SimManager sim_manager(param_file, farm_data_file);
		sim_manager.runSims(num_reps, output_file, overrides);
		
		auto time_end = time(NULL);
			
		std::cout << "\rProgress: COMPLETE - Time taken: " << (time_end - time_start)/3600.0 << " hours" << std::endl;
	}

}
