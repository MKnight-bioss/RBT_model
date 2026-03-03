#ifndef FARM_H
#define FARM_H

#include <vector>

struct Farm{
	int id;
	
	int disease_state = 0;
	bool first_inf = false;
	bool recovered = false;
	double inf_time = -1.0;
	
	double eta_equi;
	double zeta_equi;
	double eta;
	double zeta;
	
	double demand = 0.0;
	double supply = 0.0;
	double fractional_demand = 0.0;
	double fractional_supply = 0.0;
	double floor_demand = 0.0;
	double floor_supply = 0.0;
	
	double b;
	double a;
	double d;
	
	std::vector<int> sellers;
	std::vector<int> buyers;
	std::vector<int> distinct_sellers;
	std::vector<double> seller_supply_pow_m;
	
	int position_in_rate_vector;
	double total_formation_rate = 0.0;
	double total_cessation_rate = 0.0;
	double total_trade_rate = 0.0;
	
	double risk_score = 0.0;
	bool is_testing = false;
	bool is_sharing = true;
	bool to_be_set_max_risk_score = false;
	
	double sum_demand_equi = 0.0;
	double sum_supply_equi = 0.0;
	double sum_sellers_equi = 0.0;
	double sum_distinct_sellers_equi = 0.0;
	double sum_num_trades_equi = 0.0;
	double sum_batch_size_equi = 0.0;
	double sum_in_flow_equi = 0.0;
	double sum_out_flow_equi = 0.0;
	double sum_batches_rejected_equi = 0.0;
	double sum_risk_score_equi = 0.0;
	double sum_num_sales_equi = 0.0;
	
	double sum_demand_equi_change = 0.0;
	double sum_supply_equi_change = 0.0;
	double sum_sellers_equi_change = 0.0;
	double sum_distinct_sellers_equi_change = 0.0;
	double sum_risk_score_equi_change = 0.0;
	
	double num_times_inf_equi = 0.0;
	double tot_num_inf_trades_equi = 0.0;
	double tot_inf_animals_equi = 0.0;
	double tot_missed_animals_equi = 0.0;
	double tot_time_inf_equi = 0.0;
};

#endif