#include "SimManager.h"
#include "Sim.h"

#include <iostream>
#include <fstream>
#include <omp.h>
#include <vector>

struct FarmEquis{
	
	int id;
	double demand_rate;
	double supply_rate;
	double demand;
	double supply;
	double sellers;
	double distinct_sellers;
	double num_trades;
	double batch_size;
	double in_flow;
	double out_flow;
	double batches_rejected;
	double risk_score;
	double times_inf;
	double num_inf_trades;
	double tot_inf_animals;
	double tot_missed_animals;
	double tot_time_inf;
	double num_sales;
};

struct SimResult{
	
	int vector_size;
	double increment_size;
	
	double N;
	double max_risk_score;
	double test_sensitivity;
	double frac_testing_purchased_batches;
	double frac_not_sharing;
	double risk_score_not_sharing;
	double equi_begin_time;
	double equi_end_time;
	
	std::vector<double> demand_rate_unit_time;
	std::vector<double> supply_rate_unit_time;
	std::vector<double> demand_unit_time;
	std::vector<double> supply_unit_time;
	std::vector<double> sellers_unit_time;
	std::vector<double> distinct_sellers_unit_time;
	std::vector<double> num_trades_unit_time;
	std::vector<double> batch_size_unit_time;
	std::vector<double> in_vol_unit_time;
	std::vector<double> price_unit_time;
	std::vector<double> risk_score_unit_time;
	std::vector<double> infected_unit_time;
	std::vector<double> rejected_batches_unit_time;
	std::vector<double> animals_removed_unit_time;
	std::vector<double> inf_animals_unit_time;
	std::vector<double> detected_animals_unit_time;
	std::vector<double> missed_animals_unit_time;
	std::vector<double> pct_animals_removed_unit_time;
	std::vector<double> pct_inf_animals_unit_time;
	std::vector<double> pct_detected_animals_unit_time;
	std::vector<double> pct_missed_animals_unit_time;
	std::vector<double> inf_trades_unit_time;
	
	std::vector<FarmEquis> farm_equis;
};

SimResult SimManager::runSingleSim(int sim_id, const std::unordered_map<std::string, double>& overrides){
	
	bool sim_succeeded = false;
	
	SimResult result;
	
	while (!sim_succeeded){
		
		Sim sim;
		sim.initialiseRNG(sim_id);
		sim.readSimParams(param_file);
		sim.initialiseFarms(farm_data_file);
		sim.overrideParams(overrides);
		sim.setVectorSizes();
		sim.setBlockVectors();
		sim.setInitRates();
		
		try{
			sim.runSim();
			sim_succeeded = true;
		} catch (const std::runtime_error& e){
			
			continue;
		}
		
		sim.formatData(false);
	
		result.vector_size = sim.vector_size;
		result.increment_size = sim.increment_size;
		result.N = sim.N;
		result.max_risk_score = sim.max_risk_score;
		result.test_sensitivity = sim.test_sensitivity;
		result.frac_testing_purchased_batches = sim.frac_testing_purchased_batches;
		result.frac_not_sharing = sim.frac_not_sharing;
		result.risk_score_not_sharing = sim.risk_score_not_sharing;
		result.equi_begin_time = sim.equi_begin_time;
		result.equi_end_time = sim.equi_end_time;
		
		result.demand_rate_unit_time = sim.demand_rate_unit_time;
		result.supply_rate_unit_time = sim.supply_rate_unit_time;
		result.demand_unit_time = sim.demand_unit_time;
		result.supply_unit_time = sim.supply_unit_time;
		result.sellers_unit_time = sim.sellers_unit_time;
		result.distinct_sellers_unit_time = sim.distinct_sellers_unit_time;
		result.num_trades_unit_time = sim.num_trades_unit_time;
		result.batch_size_unit_time = sim.batch_size_unit_time;
		result.in_vol_unit_time = sim.in_vol_unit_time;
		result.price_unit_time = sim.price_unit_time;
		result.risk_score_unit_time = sim.risk_score_unit_time;
		result.infected_unit_time = sim.infected_unit_time;
		result.rejected_batches_unit_time = sim.rejected_batches_unit_time;
		result.animals_removed_unit_time = sim.animals_removed_unit_time;
		result.inf_animals_unit_time = sim.inf_animals_unit_time;
		result.detected_animals_unit_time = sim.detected_animals_unit_time;
		result.missed_animals_unit_time = sim.missed_animals_unit_time;
		result.pct_animals_removed_unit_time = sim.pct_animals_removed_unit_time;
		result.pct_inf_animals_unit_time = sim.pct_inf_animals_unit_time;
		result.pct_detected_animals_unit_time = sim.pct_detected_animals_unit_time;
		result.pct_missed_animals_unit_time = sim.pct_missed_animals_unit_time;
		result.inf_trades_unit_time = sim.inf_trades_unit_time;
		
		result.farm_equis.reserve(sim.farms.size());
		
		for (const auto& farm : sim.farms){
			
			FarmEquis equi;
			equi.id = farm.id;
			equi.demand_rate = farm.eta_equi;
			equi.supply_rate = farm.zeta_equi;
			equi.demand = farm.sum_demand_equi;
			equi.supply = farm.sum_supply_equi;
			equi.sellers = farm.sum_sellers_equi;
			equi.distinct_sellers = farm.sum_distinct_sellers_equi / (sim.equi_end_time - sim.equi_begin_time);
			equi.num_trades = farm.sum_num_trades_equi / (sim.equi_end_time - sim.equi_begin_time);
			equi.batch_size = farm.sum_batch_size_equi;
			equi.in_flow = farm.sum_in_flow_equi / (sim.equi_end_time - sim.equi_begin_time);
			equi.out_flow = farm.sum_out_flow_equi / (sim.equi_end_time - sim.equi_begin_time);
			equi.batches_rejected = farm.sum_batches_rejected_equi;
			equi.risk_score = farm.sum_risk_score_equi;
			equi.times_inf = farm.num_times_inf_equi;
			equi.num_inf_trades = farm.tot_num_inf_trades_equi / (sim.equi_end_time - sim.equi_begin_time);
			equi.tot_inf_animals = farm.tot_inf_animals_equi / (sim.equi_end_time - sim.equi_begin_time);
			equi.tot_missed_animals = farm.tot_missed_animals_equi / (sim.equi_end_time - sim.equi_begin_time);
			equi.tot_time_inf = farm.tot_time_inf_equi / (sim.equi_end_time - sim.equi_begin_time);
			equi.num_sales = farm.sum_num_sales_equi / (sim.equi_end_time - sim.equi_begin_time);
			result.farm_equis.push_back(equi);
		}
		
	}
	
	return result;
}

SimManager::SimManager(const std::string& paramFile, const std::string& farmDataFile) : param_file(paramFile), farm_data_file(farmDataFile) {}

void SimManager::runSims(int numSims, const std::string& outputFile, const std::unordered_map<std::string, double>& overrides){
	
	std::vector<SimResult> results(numSims);
	
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < numSims; ++i){
		
		results[i] = runSingleSim(i, overrides);
	}
	
	std::ofstream output_file(outputFile);
	if (!output_file.is_open()){
		throw std::runtime_error("Could not open output file");
	}
	
	int num_farms = results[0].N;
	std::vector<double> demand_rate(num_farms, 0.0);
	std::vector<double> supply_rate(num_farms, 0.0);
	std::vector<double> demand(num_farms, 0.0);
	std::vector<double> supply(num_farms, 0.0);
	std::vector<double> sellers(num_farms, 0.0);
	std::vector<double> distinct_sellers(num_farms, 0.0);
	std::vector<double> num_trades(num_farms, 0.0);
	std::vector<double> batch_size(num_farms, 0.0);
	std::vector<double> in_flow(num_farms, 0.0);
	std::vector<double> out_flow(num_farms, 0.0);
	std::vector<double> batches_rejected(num_farms, 0.0);
	std::vector<double> risk_score(num_farms, 0.0);
	std::vector<double> num_times_inf(num_farms, 0.0);
	std::vector<double> tot_num_inf_trades(num_farms, 0.0);
	std::vector<double> tot_inf_animals(num_farms, 0.0);
	std::vector<double> tot_missed_animals(num_farms, 0.0);
	std::vector<double> tot_time_inf(num_farms, 0.0);
	std::vector<double> num_sales(num_farms, 0.0);
	
	for (const auto& res : results){
		for (int i = 0; i < num_farms; ++i){
			demand_rate[i] = res.farm_equis[i].demand_rate;
			supply_rate[i] = res.farm_equis[i].supply_rate;
			demand[i] += res.farm_equis[i].demand / numSims;
			supply[i] += res.farm_equis[i].supply / numSims;
			sellers[i] += res.farm_equis[i].sellers / numSims;
			distinct_sellers[i] += res.farm_equis[i].distinct_sellers / numSims;
			num_trades[i] += res.farm_equis[i].num_trades / numSims;
			batch_size[i] += res.farm_equis[i].batch_size / numSims;
			in_flow[i] += res.farm_equis[i].in_flow / numSims;
			out_flow[i] += res.farm_equis[i].out_flow / numSims;
			batches_rejected[i] += res.farm_equis[i].batches_rejected / numSims;
			risk_score[i] += res.farm_equis[i].risk_score / numSims;
			num_times_inf[i] += res.farm_equis[i].times_inf / numSims;
			tot_num_inf_trades[i] += res.farm_equis[i].num_inf_trades / numSims;
			tot_inf_animals[i] += res.farm_equis[i].tot_inf_animals / numSims;
			tot_missed_animals[i] += res.farm_equis[i].tot_missed_animals / numSims;
			tot_time_inf[i] += res.farm_equis[i].tot_time_inf / numSims;
			num_sales[i] += res.farm_equis[i].num_sales / numSims;
		}
	}
	
	/*output_file << "farm,demand_rate,supply_rate,demand,supply,sellers,distinct_sellers,num_trades,batch_size,in_flow,out_flow,batches_rejected,risk_score,num_times_inf,num_inf_trades,inf_animals,missed_animals,time_inf,sales" << std::endl;
	
	for (int i = 0; i < num_farms; ++i){
		output_file << results[0].farm_equis[i].id << "," << demand_rate[i] << "," << supply_rate[i] << "," << demand[i] << "," << supply[i] << "," << sellers[i] << "," 
					<< distinct_sellers[i] << "," << num_trades[i] << "," << batch_size[i] << "," << in_flow[i] << "," << out_flow[i] << "," << batches_rejected[i] << "," 
					<< risk_score[i] << "," << num_times_inf[i] << "," << tot_num_inf_trades[i] << "," << tot_inf_animals[i] << "," << tot_missed_animals[i] << "," << tot_time_inf[i] << "," << num_sales[i] << std::endl;
	}*/
	
	output_file << "sim,time,max_risk_score,test_sensitivity,frac_testing,frac_not_sharing,risk_score_not_sharing,eta,zeta,price,demand,supply,sellers,distinct_sellers,num_trades,batch_size,in_flow,infected,risk_score,rejected_batches,animals_removed,inf_animals,detected_animals,missed_animals,pct_animals_removed,pct_inf_animals,pct_detected_animals,pct_missed_animals,inf_trades" << std::endl;
	
	for (int i = 0; i < numSims; ++i){
		auto &result = results[i];
		
		for (int j = 0; j < result.vector_size; ++j){
			
			output_file << i + 1 << "," << j * result.increment_size << "," << result.max_risk_score << "," << result.test_sensitivity << "," << result.frac_testing_purchased_batches << "," << result.frac_not_sharing << "," << result.risk_score_not_sharing << "," << result.demand_rate_unit_time[j] << "," << result.supply_rate_unit_time[j] << 
					"," << result.price_unit_time[j] << "," << result.demand_unit_time[j] << "," << result.supply_unit_time[j] << "," << result.sellers_unit_time[j] << "," << result.distinct_sellers_unit_time[j] << "," << result.num_trades_unit_time[j]/(result.increment_size * result.N) << 
					"," << result.batch_size_unit_time[j] << "," << result.in_vol_unit_time[j]/(result.increment_size * result.N) << "," << result.infected_unit_time[j] << "," << result.risk_score_unit_time[j] << "," << result.rejected_batches_unit_time[j]/(result.increment_size * result.N) <<
					"," << result.animals_removed_unit_time[j] << "," << result.inf_animals_unit_time[j] << "," << result.detected_animals_unit_time[j] << "," << result.missed_animals_unit_time[j] << "," << result.pct_animals_removed_unit_time[j] << "," << result.pct_inf_animals_unit_time[j] << "," << result.pct_detected_animals_unit_time[j] << "," << result.pct_missed_animals_unit_time[j] <<
					"," << result.inf_trades_unit_time[j] << std::endl;
		}
	}
}