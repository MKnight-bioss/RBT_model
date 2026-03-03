#include "Sim.h"
#include <cmath>
#include <algorithm>
#include <vector> 
#include <sstream>

void Sim::initialiseRNG(int simID){
	rng = gsl_rng_alloc(gsl_rng_mt19937);
	unsigned long int seed = time(NULL) + simID;
	gsl_rng_set(rng, seed);
}

void Sim::createOutputFile(const std::string& file_name){
	output_file.open(file_name);
	
	if (!output_file.is_open()){
		std::cerr << "File not found" << std::endl;
		return;
	}
	
	output_file << "time,max_risk_score,test_sensitivity,frac_testing,eta,zeta,price,demand,supply,sellers,distinct_sellers,num_trades,batch_size,in_flow,infected,risk_score,rejected_batches,animals_removed,inf_animals,detected_animals,missed_animals,pct_animals_removed,pct_inf_animals,pct_detected_animals,pct_missed_animals,inf_trades" << std::endl;
	output_file << 0.0 << "," << max_risk_score << "," << test_sensitivity << "," << frac_testing_purchased_batches << "," << avg_eta_init << "," << avg_zeta_init << "," << price << "," 
		<< 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << "," << 0.0 << std::endl;
}

void Sim::formatData(bool write_to_file){
	demand_rate_unit_time[0] = avg_eta_init;
	supply_rate_unit_time[0] = avg_zeta_init;
	price_unit_time[0] = 1.0;
	
	transform(in_vol_unit_time.begin(), in_vol_unit_time.end(), num_trades_unit_time.begin(), batch_size_unit_time.begin(), std::divides<double>());
	
	for (int i = 0; i < N; ++i){
		if (farms[i].sum_demand_equi_change > 0.0)
			farms[i].sum_demand_equi /= farms[i].sum_demand_equi_change;
		
		if (farms[i].sum_supply_equi_change > 0.0)
			farms[i].sum_supply_equi /= farms[i].sum_supply_equi_change;
		
		if (farms[i].sum_sellers_equi_change > 0.0)
			farms[i].sum_sellers_equi /= farms[i].sum_sellers_equi_change;
		
		if (farms[i].sum_num_trades_equi > 0.0)
			farms[i].sum_batch_size_equi /= farms[i].sum_num_trades_equi;
		
		if (farms[i].sum_distinct_sellers_equi_change > 0.0)
			farms[i].sum_distinct_sellers_equi /= farms[i].sum_distinct_sellers_equi_change;
		
		if (farms[i].sum_num_trades_equi > 0.0)
			farms[i].sum_batches_rejected_equi /= farms[i].sum_num_trades_equi;
		
		if (farms[i].sum_risk_score_equi_change > 0.0)
			farms[i].sum_risk_score_equi /= farms[i].sum_risk_score_equi_change;
	}
	
	for (int i = 0; i < vector_size; ++i){
		if (i > 0){
			demand_rate_unit_time[i] += demand_rate_unit_time[i-1];
			supply_rate_unit_time[i] += supply_rate_unit_time[i-1];
			demand_unit_time[i] += demand_unit_time[i-1];
			supply_unit_time[i] += supply_unit_time[i-1];
			sellers_unit_time[i] += sellers_unit_time[i-1];
			if (num_trades_unit_time[i] == 0.0)
				batch_size_unit_time[i] = batch_size_unit_time[i-1];
			price_unit_time[i] += price_unit_time[i-1];
			risk_score_unit_time[i] += risk_score_unit_time[i-1];
			infected_unit_time[i] += infected_unit_time[i-1];
			
			if (num_trades_unit_time[i] > 0.0){
				animals_removed_unit_time[i] /= num_trades_unit_time[i];
				inf_animals_unit_time[i] /= num_trades_unit_time[i];
				detected_animals_unit_time[i] /= num_trades_unit_time[i];
				missed_animals_unit_time[i] /= num_trades_unit_time[i];
			}
			
			if (inf_trades_unit_time[i] > 0.0){
				pct_animals_removed_unit_time[i] /= inf_trades_unit_time[i];
				pct_inf_animals_unit_time[i] /= inf_trades_unit_time[i];
				pct_detected_animals_unit_time[i] /= inf_trades_unit_time[i];
				pct_missed_animals_unit_time[i] /= inf_trades_unit_time[i];
			}
		}
		if (write_to_file){
			
			output_file << i * increment_size << "," << max_risk_score << "," << test_sensitivity << "," << frac_testing_purchased_batches << "," << demand_rate_unit_time[i] << "," << supply_rate_unit_time[i] << 
				"," << price_unit_time[i] << "," << demand_unit_time[i] << "," << supply_unit_time[i] << "," << sellers_unit_time[i] << "," << distinct_sellers_unit_time[i] << "," << num_trades_unit_time[i] / increment_size << 
				"," << batch_size_unit_time[i] << "," << in_vol_unit_time[i] / increment_size << "," << infected_unit_time[i] << "," << risk_score_unit_time[i] << "," << rejected_batches_unit_time[i] / increment_size << std::endl;
		}
			
	}
	
	output_file.close();
}

void Sim::readSimParams(const std::string& file_name){
	
	std::ifstream file(file_name);
	
	if (!file.is_open()){
		std::cerr << "Sim param file not found" << std::endl;
		return;
	}
	
	std::string line;
	std::getline(file, line);
	
	while (std::getline(file, line)){
		std::istringstream iss(line);
		
		std::vector<std::string> tokens;
		std::string token;
		
		while (getline(iss, token, ','))
			tokens.push_back(token);
		
		if ((std::string) tokens[0] == "per_animal_inf_prob") per_animal_inf_prob = stod(tokens[1]);
		else if ((std::string) tokens[0] == "gamma") gamma = stod(tokens[1]);
		else if ((std::string) tokens[0] == "price_volatility") price_volatility = stod(tokens[1]);
		else if ((std::string) tokens[0] == "demand_elasticity") demand_elasticity = stod(tokens[1]);
		else if ((std::string) tokens[0] == "supply_elasticity") supply_elasticity = stod(tokens[1]);
		else if ((std::string) tokens[0] == "t_max") t_max = stod(tokens[1]);
		else if ((std::string) tokens[0] == "disease_intro_time") disease_intro_time = stod(tokens[1]);
		else if ((std::string) tokens[0] == "control_intro_time") control_intro_time = stod(tokens[1]);
		else if ((std::string) tokens[0] == "frac_testing_purchased_batches") frac_testing_purchased_batches = stod(tokens[1]);
		else if ((std::string) tokens[0] == "test_sensitivity") test_sensitivity = stod(tokens[1]);
		else if ((std::string) tokens[0] == "max_risk_score") max_risk_score = stod(tokens[1]);
		else if ((std::string) tokens[0] == "time_to_min_risk_score") time_to_min_risk_score = stod(tokens[1]);
		else if ((std::string) tokens[0] == "frac_not_sharing") frac_not_sharing = stod(tokens[1]);
		else if ((std::string) tokens[0] == "risk_score_not_sharing") risk_score_not_sharing = stod(tokens[1]);
		else if ((std::string) tokens[0] == "equi_begin_time") equi_begin_time = stod(tokens[1]);
		else if ((std::string) tokens[0] == "equi_end_time") equi_end_time = stod(tokens[1]);
	}
	file.close();
}

void Sim::initialiseFarms(const std::string& file_name){
	std::ifstream file(file_name);
	
	if (!file.is_open()){
		std::cerr << "Farm param file not found" << std::endl;
		return;
	}
	
	std::string line;
	std::getline(file,line);
	
	while (std::getline(file, line)){
		std::istringstream iss(line);
		
		std::vector<std::string> tokens;
		std::string token;
		
		while (getline(iss, token, ','))
			tokens.push_back(token);
		
		Farm f;
		f.id = stoi(tokens[0]);
		f.eta_equi = stod(tokens[1]);
		f.eta = stod(tokens[1]);
		f.zeta_equi = stod(tokens[2]);
		f.zeta = stod(tokens[2]);
		f.d = 1.0 / stod(tokens[4]);
		f.a = stod(tokens[6]);
		f.b = stod(tokens[7]);
		farms.push_back(f);
	}
	
	N = static_cast<double>(farms.size());
	S = N;
	I = N - S;
}

void Sim::overrideParams(const std::unordered_map<std::string, double>& overrides){
	for (auto& [param_name, val] : overrides){
		if (param_name == "per_animal_inf_prob") per_animal_inf_prob = val;
		else if (param_name == "gamma") gamma = val;
		else if (param_name == "price_volatility") price_volatility = val;
		else if (param_name == "demand_elasticity") demand_elasticity = val;
		else if (param_name == "supply_elasticity") supply_elasticity = val;
		else if (param_name == "disease_intro_time") disease_intro_time = val;
		else if (param_name == "control_intro_time") control_intro_time = val;
		else if (param_name == "frac_testing_purchased_batches") frac_testing_purchased_batches = val;
		else if (param_name == "test_sensitivity") test_sensitivity = val;
		else if (param_name == "max_risk_score") max_risk_score = val;
		else if (param_name == "time_to_min_risk_score") time_to_min_risk_score = val;
		else if (param_name == "frac_not_sharing") frac_not_sharing = val;
		else if (param_name == "risk_score_not_sharing") risk_score_not_sharing = val;
	}
}

void Sim::setVectorSizes(){
	vector_size = int(t_max / increment_size) + 1;
	
	demand_rate_unit_time.resize(vector_size, 0.0);
	supply_rate_unit_time.resize(vector_size, 0.0);
	demand_unit_time.resize(vector_size, 0.0);
	supply_unit_time.resize(vector_size, 0.0);
	sellers_unit_time.resize(vector_size, 0.0);
	distinct_sellers_unit_time.resize(vector_size, 0.0);
	num_trades_unit_time.resize(vector_size, 0.0);
	batch_size_unit_time.resize(vector_size, 0.0);
	in_vol_unit_time.resize(vector_size, 0.0);
	price_unit_time.resize(vector_size, 0.0);
	risk_score_unit_time.resize(vector_size, 0.0);
	infected_unit_time.resize(vector_size, 0.0);
	rejected_batches_unit_time.resize(vector_size, 0.0);
	animals_removed_unit_time.resize(vector_size, 0.0);
	inf_animals_unit_time.resize(vector_size, 0.0);
	detected_animals_unit_time.resize(vector_size, 0.0);
	missed_animals_unit_time.resize(vector_size, 0.0);
	pct_animals_removed_unit_time.resize(vector_size, 0.0);
	pct_inf_animals_unit_time.resize(vector_size, 0.0);
	pct_detected_animals_unit_time.resize(vector_size, 0.0);
	pct_missed_animals_unit_time.resize(vector_size, 0.0);
	inf_trades_unit_time.resize(vector_size, 0.0);
	
	farm_supply_last_update.resize((int) N, 0.0);
	farm_demand_last_update.resize((int) N, 0.0);
}

void Sim::resetVectors(){
	std::fill(demand_rate_unit_time.begin(), demand_rate_unit_time.end(), 0.0);
	std::fill(supply_rate_unit_time.begin(), supply_rate_unit_time.end(), 0.0);
	std::fill(demand_unit_time.begin(), demand_unit_time.end(), 0.0);
	std::fill(supply_unit_time.begin(), supply_unit_time.end(), 0.0);
	std::fill(sellers_unit_time.begin(), sellers_unit_time.end(), 0.0);
	std::fill(distinct_sellers_unit_time.begin(), distinct_sellers_unit_time.end(), 0.0);
	std::fill(num_trades_unit_time.begin(), num_trades_unit_time.end(), 0.0);
	std::fill(batch_size_unit_time.begin(), batch_size_unit_time.end(), 0.0);
	std::fill(in_vol_unit_time.begin(), in_vol_unit_time.end(), 0.0);
	std::fill(price_unit_time.begin(), price_unit_time.end(), 0.0);
	std::fill(risk_score_unit_time.begin(), risk_score_unit_time.end(), 0.0);
	std::fill(infected_unit_time.begin(), infected_unit_time.end(), 0.0);
	std::fill(rejected_batches_unit_time.begin(), rejected_batches_unit_time.end(), 0.0);
	std::fill(animals_removed_unit_time.begin(), animals_removed_unit_time.end(), 0.0);
	std::fill(inf_animals_unit_time.begin(), inf_animals_unit_time.end(), 0.0);
	std::fill(detected_animals_unit_time.begin(), detected_animals_unit_time.end(), 0.0);
	std::fill(missed_animals_unit_time.begin(), missed_animals_unit_time.end(), 0.0);
	std::fill(pct_animals_removed_unit_time.begin(), pct_animals_removed_unit_time.end(), 0.0);
	std::fill(pct_inf_animals_unit_time.begin(), pct_inf_animals_unit_time.end(), 0.0);
	std::fill(pct_detected_animals_unit_time.begin(), pct_detected_animals_unit_time.end(), 0.0);
	std::fill(pct_missed_animals_unit_time.begin(), pct_missed_animals_unit_time.end(), 0.0);
	std::fill(inf_trades_unit_time.begin(), inf_trades_unit_time.end(), 0.0);
	
	std::fill(farm_supply_last_update.begin(), farm_supply_last_update.end(), 0.0);
	std::fill(farm_demand_last_update.begin(), farm_demand_last_update.end(), 0.0);
}

void Sim::setBlockVectors(){
	num_partitions = (int) std::floor(std::sqrt(N));
	
	formation_rates.resize(num_partitions, 0.0);
	trade_rates.resize(num_partitions, 0.0);
	supply_pow_m_block.resize(num_partitions, 0.0);
	final_farm.resize(num_partitions, 0.0);
}

void Sim::setInitRates(){
	
	int position;
	for (int i = 0; i < (int) N; ++i){
		
		position = (int) std::floor(i / (N / num_partitions));
		farms[i].position_in_rate_vector = position;
		final_farm[position] = i;
		
		farms[i].seller_supply_pow_m.resize(num_partitions, 0.0);
		
		total_demand_rate += farms[i].eta_equi;
		total_supply_rate += farms[i].zeta_equi;
		
		if (farms[i].eta_equi > 0.0)
			farms[i].fractional_demand = gsl_rng_uniform(rng);
		if (farms[i].zeta_equi > 0.0)
			farms[i].fractional_supply = gsl_rng_uniform(rng);
	}
	
	avg_eta_init = total_demand_rate / N;
	avg_zeta_init = total_supply_rate / N;
}

void Sim::resetChangeVariables(){
	
	change_to_demand = 0.0;
	change_to_supply = 0.0;
	change_to_sellers = 0.0;
	change_to_formation_rate = 0.0;
	change_to_trade_rate = 0.0;
	change_to_demand_rate = 0.0;
	change_to_supply_rate = 0.0;
	change_to_risk_score = 0.0;
	
}

double Sim::generateTimeStep(const double& total_rate, const double& base_time_step) const{
	return total_event_rate > 0.0 ? -std::log(gsl_rng_uniform_pos(rng)) / total_event_rate : base_time_step;
}

void Sim::updateTime(const double& step_amount){
	t_curr += step_amount;
}

double Sim::generateNextEventRate(const double& total_event_rate) const{
	return gsl_rng_uniform(rng) * total_event_rate;
}

Sim::Event Sim::generateNextEvent(const double& next_event_rate) const{
	if (next_event_rate < total_formation_rate)
		return Event::PARTNER_FORMATION;
	else if (next_event_rate >= total_formation_rate && next_event_rate < total_formation_rate + total_trade_rate)
		return Event::TRADE;
	else if (next_event_rate >= total_formation_rate + total_trade_rate && next_event_rate < total_formation_rate + total_trade_rate + total_rec_rate)
		return Event::RECOVERY;
	else{
		std::cout << "Invalid event " << next_event_rate << '\t' << total_event_rate << std::endl;
		return Event::INVALID_EVENT;
	}
}

int Sim::findEventBlock(std::vector<double>& rates_vec, const double& next_event_rate, double& rate_accumulator){
	
	for (std::size_t i = 0; i < rates_vec.size(); ++i){
		rate_accumulator += rates_vec[i];
		
		if (rate_accumulator >= next_event_rate){
			
			rate_accumulator -= rates_vec[i];
			return i;
		}
	}
	
	return -1;
}

int Sim::findEventPerformer(const Event& next_event, const double& next_event_rate, double& rate_accumulator){
	
	int block_index;
	int farm_iter = 0;
	
	switch(next_event){
		
		case Event::PARTNER_FORMATION: {
			
			block_index = findEventBlock(formation_rates, next_event_rate, rate_accumulator);
			
			if (block_index > 0)
				farm_iter = final_farm[block_index - 1] + 1;
			
			for (int it = farm_iter; it < final_farm[block_index] + 1; ++it){
				
				rate_accumulator += farms[it].total_formation_rate;
				
				if (rate_accumulator >= next_event_rate){
					rate_accumulator -= farms[it].total_formation_rate;
					return it;
				}
			}
			
			break;
		}
		
		case Event::TRADE: {
			
			rate_accumulator += total_formation_rate;
			
			block_index = findEventBlock(trade_rates, next_event_rate, rate_accumulator);
			
			if (block_index > 0)
				farm_iter = final_farm[block_index - 1] + 1;
			
			for (int it = farm_iter; it < final_farm[block_index] + 1; ++it){
				
				rate_accumulator += farms[it].total_trade_rate;
				
				if (rate_accumulator >= next_event_rate){
					rate_accumulator -= farms[it].total_trade_rate;
					return it;
				}
			}
			
			break;
		}
		
		case Event::RECOVERY: {
			
			rate_accumulator += total_formation_rate + total_trade_rate;
			
			for (int i = 0; i < N; ++i){
				
				if (farms[i].disease_state == 1){
					
					rate_accumulator += calcRecRate();
					
					if (rate_accumulator >= next_event_rate)
						return i;
				}
			}
			
			break;
		}
		
		case Event::INVALID_EVENT: {
			
			return -1;
			break;
		}
	}
	
	return -1;
}

int Sim::findEventPartner(const Event& next_event, const int& event_performer, const double& next_event_rate, double& rate_accumulator){
	
	switch(next_event){
		
		case Event::PARTNER_FORMATION: {
			
			double buyer_supply_pow_m = 0.0;
			
			for (int i = 0; i < num_partitions; ++i){
				
				buyer_supply_pow_m = 0.0;
				
				if (farms[event_performer].position_in_rate_vector == i)
					buyer_supply_pow_m = (1.0 - farms[event_performer].risk_score) * pow(farm_supply_last_update[event_performer], m);
				
				rate_accumulator += farms[event_performer].a * farm_demand_last_update[event_performer] * (supply_pow_m_block[i] - farms[event_performer].seller_supply_pow_m[i] - buyer_supply_pow_m) / (N-1);
				
				if (rate_accumulator >= next_event_rate){
					
					rate_accumulator -= farms[event_performer].a * farm_demand_last_update[event_performer] * (supply_pow_m_block[i] - farms[event_performer].seller_supply_pow_m[i] - buyer_supply_pow_m) / (N-1);
					
					int block_iterator = 0;
					
					if (i > 0)
						block_iterator = final_farm[i-1] + 1;
					
					for (int j = block_iterator; j < final_farm[i] + 1; ++j){
						
						if (event_performer != j && farm_supply_last_update[j] > 0.0){
							
							if( !isSeller(event_performer, j) ){
								
								rate_accumulator += calcFormationRate(event_performer, j);
								
								if (rate_accumulator >= next_event_rate)
									return j;
							}
						}
					}
				}
			}
			
			break;
		}
		
		case Event::TRADE: {
			
			for (std::size_t j = 0; j < farms[event_performer].sellers.size(); ++j){
				
				int seller_index = farms[event_performer].sellers[j];
				
				if (event_performer != seller_index && farms[seller_index].supply > 0.0){
					
					rate_accumulator += calcTradeRate(event_performer, seller_index);
					
					if (rate_accumulator >= next_event_rate){
						return seller_index;
					}
				}
			}
			
			break;
		}
	}
	
	return -1;
}

void Sim::doFormationEvent(const int& event_performer, const int& event_partner){
	
	farms[event_performer].sellers.push_back(event_partner);
	farms[event_partner].buyers.push_back(event_performer);
	farms[event_performer].seller_supply_pow_m[farms[event_partner].position_in_rate_vector] += (1.0 - farms[event_partner].risk_score)*pow(farm_supply_last_update[event_partner], m);
	
	double formation_rate = calcFormationRate(event_performer, event_partner);
	double trade_rate = calcTradeRate(event_performer, event_partner);
	
	farms[event_performer].total_formation_rate -= formation_rate;
	farms[event_performer].total_trade_rate += trade_rate;
	
	formation_rates[farms[event_performer].position_in_rate_vector] -= formation_rate;
	trade_rates[farms[event_performer].position_in_rate_vector] += trade_rate;
	
	total_formation_rate -= formation_rate;
	total_trade_rate += trade_rate;
	total_event_rate += (trade_rate - formation_rate);
	
	if (equi_intro_check && t_curr <= equi_end_time){
		
		farms[event_performer].sum_sellers_equi += farms[event_performer].sellers.size();
		++farms[event_performer].sum_sellers_equi_change;
	}
	
	updateVector(sellers_unit_time, 1.0 / N);
}

void Sim::doTradeEvent(const int& event_performer, const int& event_partner){
	
	bool new_seller = false;
	if ( !isDistinctSeller(event_performer, event_partner) ){
		new_seller = true;
		farms[event_performer].distinct_sellers.push_back(event_partner);
	}
	
	double batch_size = std::min(farms[event_performer].demand, farms[event_partner].supply);
	double animals_through = batch_size;
	double inf_animals = 0.0;
	double detected_animals = 0.0;
	bool reject_batch = false;
	
	farms[event_partner].supply -= batch_size;
	total_supply -= batch_size;
	
	if ( farms[event_partner].disease_state == 1 )
		inf_animals = calcInfAnimals(batch_size);
	
	if (!farms[event_performer].is_testing){
		
		if ( inf_animals > 0.0 && farms[event_performer].disease_state == 0 ){
			
			farms[event_performer].disease_state = 1;
			
			++I;
			--S;
			
			total_rec_rate += calcRecRate();
			total_event_rate += calcRecRate();
			
			farms[event_performer].inf_time = t_curr;
			
			if (equi_intro_check && t_curr <= equi_end_time){
		
				farms[event_performer].num_times_inf_equi += 1.0;
			}
			
			updateVector(infected_unit_time, 1.0 / N);
		}
	}
	
	else{
		
		if (inf_animals > 0.0){
			
			detected_animals = calcDetectedAnimals(inf_animals);
			
			if ( detected_animals > 0.0 ){
				
				reject_batch = true;
				
				animals_through = 0.0;
				
				if (farms[event_partner].is_sharing) 
					farms[event_partner].to_be_set_max_risk_score = true;
			}
			
			else{
				
				if (farms[event_performer].disease_state == 0){
					
					farms[event_performer].disease_state = 1;
					
					++I;
					--S;
					
					total_rec_rate += calcRecRate();
					total_event_rate += calcRecRate();
					
					farms[event_performer].inf_time = t_curr;
					
					if (equi_intro_check && t_curr <= equi_end_time){
				
						farms[event_performer].num_times_inf_equi += 1.0;
					}
					
					updateVector(infected_unit_time, 1.0/N);
				}
			}
		}
	}
	
	farms[event_performer].demand -= animals_through;
	total_demand -= animals_through;
	
	if (reject_batch){
		
		change_to_demand_rate = 0.0;
		change_to_supply_rate = 0.0;
		
		double old_price = price;
		price = calcPrice();
		
		double new_eta;
		double new_zeta;
		for (int i = 0; i < N; ++i){
			
			if (farms[i].eta_equi > 0.0){
				
				new_eta = calcDemandRate(i);
				change_to_demand_rate += (new_eta - farms[i].eta);
				farms[i].eta = new_eta;
			}
			
			if (farms[i].zeta_equi > 0.0){
				
				new_zeta = calcSupplyRate(i);
				change_to_supply_rate += (new_zeta - farms[i].zeta);
				farms[i].zeta = new_zeta;
			}
		}
		
		updateVector(demand_rate_unit_time, change_to_demand_rate / N);
		updateVector(supply_rate_unit_time, change_to_supply_rate / N);
		updateVector(price_unit_time, price - old_price);
	}
	
	change_to_trade_rate = 0.0;
	double old_trade_rate = 0.0;
	double new_trade_rate = 0.0;
	for (std::size_t seller = 0; seller < farms[event_performer].sellers.size(); ++seller){
		
		int seller_index = farms[event_performer].sellers[seller];
		
		old_trade_rate = 0.0;
		new_trade_rate = calcTradeRate(event_performer, seller_index);
		
		if (seller_index == event_partner)
			old_trade_rate = (1.0 - farms[seller_index].risk_score) * farms[event_performer].b * std::min( farms[event_performer].demand + animals_through, farms[seller_index].supply + batch_size );
		else
			old_trade_rate = (1.0 - farms[seller_index].risk_score) * farms[event_performer].b * std::min( farms[event_performer].demand + animals_through, farms[seller_index].supply );
		
		change_to_trade_rate += (new_trade_rate - old_trade_rate);
		farms[event_performer].total_trade_rate += (new_trade_rate - old_trade_rate);
		trade_rates[farms[event_performer].position_in_rate_vector] += (new_trade_rate - old_trade_rate);
	}
	
	for (std::size_t buyer = 0; buyer < farms[event_partner].buyers.size(); ++buyer){
		
		if (event_performer != farms[event_partner].buyers[buyer]){
			
			int buyer_index = farms[event_partner].buyers[buyer];
			
			old_trade_rate = (1.0 - farms[event_partner].risk_score) * farms[buyer_index].b * std::min( farms[buyer_index].demand, farms[event_partner].supply + batch_size );
			new_trade_rate = calcTradeRate(buyer_index, event_partner);
			
			change_to_trade_rate += (new_trade_rate - old_trade_rate);
			farms[buyer_index].total_trade_rate += (new_trade_rate - old_trade_rate);
			trade_rates[farms[buyer_index].position_in_rate_vector] += (new_trade_rate - old_trade_rate);
		}
	}
	
	total_trade_rate += change_to_trade_rate;
	total_event_rate += change_to_trade_rate;
	
	if (equi_intro_check && t_curr <= equi_end_time){
		farms[event_performer].sum_batch_size_equi += animals_through;
		farms[event_performer].sum_in_flow_equi += animals_through;
		farms[event_performer].sum_demand_equi += farms[event_performer].demand;
		++farms[event_performer].sum_num_trades_equi;
		++farms[event_performer].sum_demand_equi_change;
		
		if (new_seller){
			
			farms[event_performer].sum_distinct_sellers_equi += farms[event_performer].distinct_sellers.size();
			++farms[event_performer].sum_distinct_sellers_equi_change;
		}
		
		if (reject_batch)
			++farms[event_performer].sum_batches_rejected_equi;
		
		farms[event_partner].sum_supply_equi += farms[event_partner].supply;
		farms[event_partner].sum_out_flow_equi += batch_size;
		++farms[event_partner].sum_supply_equi_change;
		++farms[event_partner].sum_num_sales_equi;
		
		if (inf_animals > 0.0){
			farms[event_performer].tot_num_inf_trades_equi += 1.0;
			farms[event_performer].tot_inf_animals_equi += inf_animals;
			farms[event_performer].tot_missed_animals_equi += (inf_animals - detected_animals);
		}
	}
	
	updateVector(demand_unit_time, -animals_through / N);
	updateVector(supply_unit_time, -batch_size / N);
	updateVector(in_vol_unit_time, animals_through);
	updateVector(num_trades_unit_time, 1.0);
	if (new_seller) 
		updateVector(distinct_sellers_unit_time, 1.0 / N);
	
	if (reject_batch)
		updateVector(rejected_batches_unit_time, 1.0);
	
	updateVector(animals_removed_unit_time, batch_size - animals_through);
	updateVector(inf_animals_unit_time, inf_animals);
	updateVector(detected_animals_unit_time, detected_animals);
	updateVector(missed_animals_unit_time, inf_animals - detected_animals);
	
	updateVector(pct_animals_removed_unit_time, 100*(batch_size - animals_through)/batch_size);
	updateVector(pct_inf_animals_unit_time, 100*inf_animals/batch_size);
	if (inf_animals > 0.0){
		updateVector(pct_detected_animals_unit_time, 100*detected_animals/inf_animals);
		updateVector(pct_missed_animals_unit_time, 100*(inf_animals - detected_animals)/inf_animals);
		updateVector(inf_trades_unit_time, 1.0);
	}
}

void Sim::doRecEvent(const int& event_performer){
	
	farms[event_performer].disease_state = 0;
	
	++S;
	--I;
	
	total_rec_rate -= calcRecRate();
	total_event_rate -= calcRecRate();
	
	if (equi_intro_check && t_curr <= equi_end_time){
		
		if (farms[event_performer].inf_time < equi_begin_time)
			farms[event_performer].tot_time_inf_equi += t_curr - equi_begin_time;
		else
			farms[event_performer].tot_time_inf_equi += t_curr - farms[event_performer].inf_time;
	}
	
	updateVector(infected_unit_time, -1.0/N);
}

void Sim::doNextEvent(const Event& next_event, const double& next_event_rate){
	double partial_rate = 0.0;
	
	int event_performer = findEventPerformer(next_event, next_event_rate, partial_rate);
	
	switch(next_event){
		
		case Event::PARTNER_FORMATION:{
			
			int event_partner = findEventPartner(next_event, event_performer, next_event_rate, partial_rate);
			doFormationEvent(event_performer, event_partner);
			return;
		}
		
		case Event::TRADE:{
			
			int event_partner = findEventPartner(next_event, event_performer, next_event_rate, partial_rate);
			doTradeEvent(event_performer, event_partner);
			return;
		}
		
		case Event::RECOVERY:{
			
			doRecEvent(event_performer);
			return;
		}
		
		case Event::INVALID_EVENT:{
			
			if (total_event_rate > 0.0){
				std::cerr << "Not a valid event" << std::endl;
				t_curr = t_max + 1.0;
			}
			
			return;
		}
		
	}
}

void Sim::updateSupplyDemand(){
	
	farms_demand_updated.clear();
	farms_supply_updated.clear();
	
	change_to_demand = 0.0;
	change_to_supply = 0.0;
	
	for (int i = 0; i < N; ++i){
		
		if (farms[i].eta > 0.0){
			
			farms[i].fractional_demand += farms[i].eta * (t_curr - t_last_demand_supply_update);
			farms[i].floor_demand = std::floor(farms[i].fractional_demand);
			
			if (farms[i].floor_demand >= 1.0){
				
				farms[i].demand += farms[i].floor_demand;
				farms[i].fractional_demand -= farms[i].floor_demand;
				
				total_demand += farms[i].floor_demand;
				change_to_demand += farms[i].floor_demand;
				
				farms_demand_updated.push_back(i);
				
				if (equi_intro_check && t_curr <= equi_end_time){
					
					farms[i].sum_demand_equi += farms[i].demand;
					++farms[i].sum_demand_equi_change;
				}
			}
		}
		
		if (farms[i].zeta > 0.0){
			
			farms[i].fractional_supply += farms[i].zeta * (t_curr - t_last_demand_supply_update);
			farms[i].floor_supply = std::floor(farms[i].fractional_supply);
			
			if (farms[i].floor_supply >= 1.0){
				
				farms[i].supply += farms[i].floor_supply;
				farms[i].fractional_supply -= farms[i].floor_supply;
				
				total_supply += farms[i].floor_supply;
				change_to_supply += farms[i].floor_supply;
				
				farms_supply_updated.push_back(i);
				
				if (equi_intro_check && t_curr <= equi_end_time){
					
					farms[i].sum_supply_equi += farms[i].supply;
					++farms[i].sum_supply_equi_change;
				}
			}
		}
	}
	
	change_to_demand_rate = 0.0;
	change_to_supply_rate = 0.0;
	
	double old_price = price;
	price = calcPrice();
	
	double new_eta;
	double new_zeta;
	for (int i = 0; i < N; ++i){
		
		if (farms[i].eta_equi > 0.0){
			
			new_eta = calcDemandRate(i);
			change_to_demand_rate += (new_eta - farms[i].eta);
			farms[i].eta = new_eta;
		}
		
		if (farms[i].zeta_equi > 0.0){
			
			new_zeta = calcSupplyRate(i);
			change_to_supply_rate += (new_zeta - farms[i].zeta);
			farms[i].zeta = new_zeta;
		}
	}
	
	change_to_trade_rate = 0.0;
	double old_trade_rate = 0.0;
	double new_trade_rate = 0.0;
	for (std::size_t i = 0; i < farms_demand_updated.size(); ++i){
		
		int farm_index = farms_demand_updated[i];
		
		for (std::size_t j = 0; j < farms[farm_index].sellers.size(); ++j){
			
			int seller_index = farms[farm_index].sellers[j];
			
			if (farms[seller_index].supply == 0.0)
				continue;
				
			old_trade_rate = (1.0 - farms[seller_index].risk_score) * farms[farm_index].b * std::min( farms[farm_index].demand - farms[farm_index].floor_demand, farms[seller_index].supply - farms[seller_index].floor_supply );
			new_trade_rate = calcTradeRate(farm_index, seller_index);
			
			change_to_trade_rate += (new_trade_rate - old_trade_rate);
			farms[farm_index].total_trade_rate += (new_trade_rate - old_trade_rate);
			trade_rates[farms[farm_index].position_in_rate_vector] += (new_trade_rate - old_trade_rate);
			
		}
	}
	
	for (std::size_t i = 0; i < farms_supply_updated.size(); ++i){
		
		int farm_index = farms_supply_updated[i];
		
		for (std::size_t j = 0; j < farms[farm_index].buyers.size(); ++j){
			
			int buyer_index = farms[farm_index].buyers[j];
			
			if ( find(farms_demand_updated.begin(), farms_demand_updated.end(), buyer_index) != farms_demand_updated.end() )
				continue;
			
			if (farms[buyer_index].demand == 0.0)
				continue;
			
			old_trade_rate = (1.0 - farms[farm_index].risk_score) * farms[buyer_index].b * std::min( farms[buyer_index].demand, farms[farm_index].supply - farms[farm_index].floor_supply );
			new_trade_rate = calcTradeRate(buyer_index, farm_index);
			
			change_to_trade_rate += (new_trade_rate - old_trade_rate);
			farms[buyer_index].total_trade_rate += (new_trade_rate - old_trade_rate);
			trade_rates[farms[buyer_index].position_in_rate_vector] += (new_trade_rate - old_trade_rate);
		}
	}
	
	total_trade_rate += change_to_trade_rate;
	total_event_rate += change_to_trade_rate;
	
	updateVector(demand_unit_time, change_to_demand / N);
	updateVector(supply_unit_time, change_to_supply / N);
	updateVector(demand_rate_unit_time, change_to_demand_rate / N);
	updateVector(supply_rate_unit_time, change_to_supply_rate / N);
	updateVector(price_unit_time, price - old_price);
	
	t_last_demand_supply_update = t_curr;
}

void Sim::removeSellers(){
	
	change_to_sellers = 0.0;
	
	for (int i = 0; i < N; ++i){
		
		std::vector<int> sellers_to_be_removed;
		
		if (farms[i].d == 0.0)
			continue;
		
		if (farms[i].sellers.empty())
			continue;
		
		for (std::size_t j = 0; j < farms[i].sellers.size(); ++j){
			
			int seller_index = farms[i].sellers[j];
			
			if ( -std::log(gsl_rng_uniform_pos(rng)) / calcCessationRate(i, seller_index) < t_curr - t_last_sellers_cessation ){
				
				sellers_to_be_removed.push_back(seller_index);
				--change_to_sellers;
				
				double trade_rate = calcTradeRate(i, seller_index);
				double formation_rate = calcFormationRate(i, seller_index);
				
				farms[i].total_trade_rate -= trade_rate;
				farms[i].total_formation_rate += formation_rate;
				trade_rates[farms[i].position_in_rate_vector] -= trade_rate;
				formation_rates[farms[i].position_in_rate_vector] += formation_rate;
				total_trade_rate -= trade_rate;
				total_formation_rate += formation_rate;
				total_event_rate += (formation_rate - trade_rate);
				
				farms[i].seller_supply_pow_m[farms[seller_index].position_in_rate_vector] -= (1.0 - farms[seller_index].risk_score) * std::pow(farm_supply_last_update[seller_index], m);
				
				
				auto buyer_position_iter = find(farms[seller_index].buyers.begin(), farms[seller_index].buyers.end(), i);
				auto buyer_position_index = distance(farms[seller_index].buyers.begin(), buyer_position_iter);
				farms[seller_index].buyers.erase(farms[seller_index].buyers.begin() + buyer_position_index);
			}
		}
		
		for (std::size_t j = 0; j < sellers_to_be_removed.size(); ++j){
			
			auto seller_position_index = distance( farms[i].sellers.begin(), find(farms[i].sellers.begin(), farms[i].sellers.end(), sellers_to_be_removed[j]) );
			farms[i].sellers.erase( farms[i].sellers.begin() + seller_position_index );
		}
		
		if (equi_intro_check && t_curr <= equi_end_time){
			
			farms[i].sum_sellers_equi += farms[i].sellers.size();
			++farms[i].sum_sellers_equi_change;
		}
	}
	
	updateVector(sellers_unit_time, change_to_sellers / N);
	
	t_last_sellers_cessation = t_curr;
}

void Sim::updateRiskScores(){
	
	change_to_risk_score = 0.0;
	total_supply_pow_m = 0.0;
	
	for (int i = 0; i < N; ++i){
		
		if (!farms[i].is_sharing){
			if (farms[i].risk_score == 0.0)
				change_to_risk_score += risk_score_not_sharing;
			farms[i].risk_score = risk_score_not_sharing;
		}
		
		if (farms[i].to_be_set_max_risk_score){
			
			change_to_risk_score += max_risk_score - farms[i].risk_score;
			farms[i].risk_score = max_risk_score;
			farms[i].to_be_set_max_risk_score = false;
		}
		else{
			
			if (farms[i].risk_score > 0.0 && farms[i].is_sharing){
			
				double risk_score_update = max_risk_score * (t_curr - t_last_risk_scores_update) / time_to_min_risk_score;
				if (risk_score_update > farms[i].risk_score)
					risk_score_update = farms[i].risk_score;
				
				change_to_risk_score -= risk_score_update;
				farms[i].risk_score -= risk_score_update;
			}
		}
		
		if (farms[i].supply > 0.0)
			total_supply_pow_m += (1.0 - farms[i].risk_score) * std::pow(farms[i].supply, m);
		
		if (equi_intro_check && t_curr <= equi_end_time){
			
			farms[i].sum_risk_score_equi += farms[i].risk_score;
			++farms[i].sum_risk_score_equi_change;
		}
	}
	
	updateVector(risk_score_unit_time, change_to_risk_score / N);
	
	t_last_risk_scores_update = t_curr;
}

void Sim::resetDistinctSellers(){
	
	t_last_distinct_sellers_update = t_curr;
	
	for (auto& farm : farms){
		
		farm.distinct_sellers.clear();
	}
}

void Sim::introduceDisease(int& init_inf){
	
	disease_intro_check = true;
	
	int largest_seller_index = 0;
	
	while (I < init_inf){
		
		for (int i = 0; i < N; ++i){
			
			if (farms[i].zeta_equi > farms[largest_seller_index].zeta_equi)
				largest_seller_index = i;
		}
		
		farms[largest_seller_index].disease_state = 1;
		farms[largest_seller_index].first_inf = true;
		
		++I;
		--S;
		
		total_rec_rate = calcRecRate();
		total_event_rate += calcRecRate();
		
		farms[largest_seller_index].inf_time = t_curr;
		
		if (equi_intro_check && t_curr <= equi_end_time){
			
			farms[largest_seller_index].num_times_inf_equi += 1.0;
		}
		
		updateVector(infected_unit_time, 1.0/N);
	}
}

void Sim::introduceControl(){
	
	control_intro_check = true;
	
	if (frac_not_sharing > 0.0){
		std::vector<std::pair<int, double>> zetas;
		zetas.reserve(N);
		
		for (int i = 0; i < N; ++i){
			
			zetas.emplace_back(i, farms[i].zeta_equi);
		}
		
		std::sort(zetas.begin(), zetas.end(), [](std::pair<int, double> &a, std::pair<int, double> &b) {return a.second > b.second;});
		
		for (int i = 0; i < int(floor(N * frac_not_sharing)); ++i){
			
			farms[zetas[i].first].is_sharing = false;
		}
	}
	
	for (int i = 0; i < int(floor(N * frac_testing_purchased_batches)); ++i){
		
		farms[i].is_testing = true;
	}
}

void Sim::recalculateRates(){
	
	total_event_rate -= (total_formation_rate + total_trade_rate);
	total_formation_rate = 0.0;
	total_trade_rate = 0.0;
	
	std::fill(formation_rates.begin(), formation_rates.end(), 0.0);
	std::fill(trade_rates.begin(), trade_rates.end(), 0.0);
	std::fill(supply_pow_m_block.begin(), supply_pow_m_block.end(), 0.0);
	
	double seller_total_supply_pow_m;
	double sellers_supply_pow_m;
	double farms_own_supply_pow_m;
	for (int i = 0; i < N; ++i){
		
		seller_total_supply_pow_m = 0.0;
		
		farms[i].total_formation_rate = 0.0;
		farms[i].total_trade_rate = 0.0;
		
		std::fill(farms[i].seller_supply_pow_m.begin(), farms[i].seller_supply_pow_m.end(), 0.0);
		
		farms_own_supply_pow_m = (1.0 - farms[i].risk_score) * std::pow(farms[i].supply, m);
		supply_pow_m_block[farms[i].position_in_rate_vector] += farms_own_supply_pow_m;
		
		for (std::size_t j = 0; j < farms[i].sellers.size(); ++j){
			
			int seller_index = farms[i].sellers[j];
			
			sellers_supply_pow_m = (1.0 - farms[seller_index].risk_score) * std::pow(farms[seller_index].supply, m);
			seller_total_supply_pow_m += sellers_supply_pow_m;
			farms[i].seller_supply_pow_m[farms[seller_index].position_in_rate_vector] += sellers_supply_pow_m;
			farms[i].total_trade_rate += calcTradeRate(i, seller_index);
		}
		
		trade_rates[farms[i].position_in_rate_vector] += farms[i].total_trade_rate;
		total_trade_rate += farms[i].total_trade_rate;
		
		farms[i].total_formation_rate = farms[i].a * farms[i].demand * (total_supply_pow_m - seller_total_supply_pow_m - farms_own_supply_pow_m) / (N-1);
		formation_rates[farms[i].position_in_rate_vector] += farms[i].total_formation_rate;
		total_formation_rate += farms[i].total_formation_rate;
		
		farm_supply_last_update[i] = farms[i].supply;
		farm_demand_last_update[i] = farms[i].demand;
	}
	
	total_event_rate += total_formation_rate + total_trade_rate;
}

double Sim::calcFormationRate(const int& buyer, const int& seller){
	return (1.0 - farms[seller].risk_score) * farms[buyer].a * farm_demand_last_update[buyer] * pow(farm_supply_last_update[seller], m) / (N-1);
}

double Sim::calcCessationRate(const int& buyer, const int& seller){
	return farms[buyer].d / (1.0 - farms[seller].risk_score);
}

double Sim::calcTradeRate(const int& buyer, const int& seller){
	return (1.0 - farms[seller].risk_score) * farms[buyer].b * std::min(farms[buyer].demand, farms[seller].supply);
}

double Sim::calcRecRate(){
	return gamma;
}

double Sim::calcInfProb(const double& batch_size){
	return 1.0 - pow(1.0 - per_animal_inf_prob, batch_size);
}

double Sim::calcInfAnimals(const double& batch_size){
	return gsl_ran_binomial(rng, per_animal_inf_prob, batch_size);
}

bool Sim::canTransmit(const int& infectee, const int& infector){
	
	if (farms[infectee].disease_state == 0 && farms[infector].disease_state == 1)
		return true;
	
	return false;
}

double Sim::calcDetectedAnimals(const double& inf_animals){
	
	if (inf_animals > 0.0)
		return gsl_ran_binomial(rng, test_sensitivity, inf_animals);
	
	return 0.0;
}

double Sim::calcDemandRate(const int& farm){
	
	if (farms[farm].eta_equi > 0.0)
		return farms[farm].eta_equi * std::pow(price / equi_price, -demand_elasticity);
	
	return 0.0;
}

double Sim::calcSupplyRate(const int& farm){
	
	if (farms[farm].zeta_equi > 0.0)
		return farms[farm].zeta_equi * std::pow(price / equi_price, supply_elasticity);
	
	return 0.0;
}

double Sim::calcPrice(){
	
	return equi_price * std::exp(price_volatility * (total_demand - total_supply));
}

bool Sim::isSeller(const int& buyer, const int& seller){
	if (!farms[buyer].sellers.empty()){
		
		if (std::find(farms[buyer].sellers.begin(), farms[buyer].sellers.end(), seller) != farms[buyer].sellers.end())
			return true;
	}
	return false;
}

bool Sim::isDistinctSeller(const int& buyer, const int& seller){
	if (!farms[buyer].distinct_sellers.empty()){
		
		if (std::find(farms[buyer].distinct_sellers.begin(), farms[buyer].distinct_sellers.end(), seller) != farms[buyer].distinct_sellers.end())
			return true;
	}
	
	return false;
}

void Sim::updateVector(std::vector<double>& vec, const double& update_amount){
	auto time_index = int(std::floor(t_curr / increment_size)) + 1;
	if (time_index <= t_max / increment_size) 
		vec[time_index] += update_amount;
}

void Sim::runSim(){
	
	double t_last_progress_output = 0.0;
	
	while(t_curr < t_max){
		
		resetChangeVariables();
		
		if (t_curr >= equi_begin_time && !equi_intro_check){
			equi_intro_check = true;
			
			for (int i = 0; i < N; ++i){
				
				farms[i].sum_demand_equi += farms[i].demand;
				farms[i].sum_supply_equi += farms[i].supply;
				farms[i].sum_sellers_equi += farms[i].sellers.size();
				farms[i].sum_risk_score_equi += farms[i].risk_score;
				
				++farms[i].sum_demand_equi_change;
				++farms[i].sum_supply_equi_change;
				++farms[i].sum_sellers_equi_change;
				++farms[i].sum_risk_score_equi_change;
			}
		}
		
		if (t_curr >= equi_end_time && !equi_end_check){
			equi_end_check = true;
			
			for (int i = 0; i < N; ++i){
				if (farms[i].disease_state == 1){
					if (farms[i].inf_time < equi_begin_time)
						farms[i].tot_time_inf_equi += t_curr - equi_begin_time;
					else
						farms[i].tot_time_inf_equi += t_curr - farms[i].inf_time;
				}
			}
		}
		
		if (disease_intro_check && !control_intro_check && I == 0.0)
			throw std::runtime_error("Infection died out early");
		
		if (t_curr - t_last_distinct_sellers_update >= t_step_distinct_sellers_update)
			resetDistinctSellers();
		
		if (t_curr >= control_intro_time && !control_intro_check)
			introduceControl();
		
		if (t_curr >= disease_intro_time && !disease_intro_check)
			introduceDisease(I_init);
		
		if (t_curr - t_last_demand_supply_update >= t_step_demand_supply_update)
			updateSupplyDemand();
		
		if (t_curr - t_last_sellers_cessation >= t_step_sellers_cessation)
			removeSellers();
		
		if (t_curr - t_last_risk_scores_update >= t_step_risk_scores_update){
			
			updateRiskScores();
			recalculateRates();
		}
		
		const double time_jump = generateTimeStep(total_event_rate, t_step_no_event);
		
		if (time_jump < t_step_no_event){
			
			updateTime(time_jump);
			
			const double event_rate = generateNextEventRate(total_event_rate);
			const Sim::Event event = generateNextEvent(event_rate);
			
			doNextEvent(event, event_rate);
		}
		else
			updateTime(t_step_no_event);
	}
}