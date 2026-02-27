#ifndef SIM_H
#define SIM_H

#include "Farm.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

struct Sim{
	enum class Event { PARTNER_FORMATION, TRADE, RECOVERY, INVALID_EVENT };
	
	gsl_rng* rng;
	void initialiseRNG(int simID);
	
	std::ofstream output_file;
	void createOutputFile(const std::string& file_name);
	void formatData(bool write_to_file);
	
	void readSimParams(const std::string& file_name);
	
	double N;
	double S;
	double I = 0;
	int I_init = 1;
	
	double per_animal_inf_prob;
	double batch_inf_prob;
	double gamma;
	
	double m = 0.75; 
	
	double avg_eta_init = 0.0;
	double avg_zeta_init = 0.0;
	
	double equi_price = 1.0;
	double price = equi_price;
	double price_volatility;
	double demand_elasticity;
	double supply_elasticity;
	
	double total_demand = 0.0;
	double total_supply = 0.0;
	double total_supply_pow_m = 0.0;
	double change_to_demand = 0.0;
	double change_to_supply = 0.0;
	
	double change_to_sellers = 0.0;
	
	std::vector<int> farms_demand_updated;
	std::vector<int> farms_supply_updated;
	std::vector<double> farm_supply_last_update;
	std::vector<double> farm_demand_last_update;
	
	double total_event_rate = 0.0;
	double total_formation_rate = 0.0;
	double total_trade_rate = 0.0;
	double total_rec_rate = 0.0;
	double total_demand_rate = 0.0;
	double total_supply_rate = 0.0;
	
	double change_to_formation_rate = 0.0;
	double change_to_trade_rate = 0.0;
	double change_to_demand_rate = 0.0;
	double change_to_supply_rate = 0.0;
	
	int num_partitions;
	std::vector<double> formation_rates;
	std::vector<double> trade_rates;
	std::vector<double> supply_pow_m_block;
	std::vector<int> final_farm; //stores index of last farm in each rate block
	
	double t_max;
	double t_curr = 0.0;
	double t_last_distinct_sellers_update = 0.0;
	double t_last_demand_supply_update = 0.0;
	double t_last_sellers_cessation = 0.0;
	double t_last_risk_scores_update = 0.0;
	
	double t_step_distinct_sellers_update = 1.0;
	double t_step_demand_supply_update = 0.0002;
	double t_step_sellers_cessation = 0.01;
	double t_step_risk_scores_update = 0.01;
	double t_step_no_event = t_step_demand_supply_update;
	
	double disease_intro_time;
	double control_intro_time;
	double equi_begin_time;
	double equi_end_time;
	bool disease_intro_check = false;
	bool control_intro_check = false;
	bool equi_intro_check = false;
	bool equi_end_check = false;
	
	double frac_testing_purchased_batches;
	double test_sensitivity;
	double max_risk_score;
	double time_to_min_risk_score;
	double change_to_risk_score = 0.0;
	double frac_not_sharing;
	double risk_score_not_sharing;
	
	double increment_size = 0.1;
	int vector_size;
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
	
	std::vector<Farm> farms;
	void initialiseFarms(const std::string& file_name);
	void overrideParams(const std::unordered_map<std::string, double>& overrides);
	
	void setVectorSizes();
	void resetVectors();
	
	void setBlockVectors();
	void resetBlockVectors();
	
	void setInitRates();
	
	void resetChangeVariables();
	
	double generateTimeStep(const double& total_rate, const double& base_time_step) const;
	void updateTime(const double& step_amount);
	
	double generateNextEventRate(const double& total_rate) const;
	Event generateNextEvent(const double& next_event_rate) const;
	int findEventBlock(std::vector<double>& rates_vec, const double& next_event_rate, double& rate_accumulator);
	int findEventPerformer(const Event& next_event, const double& next_event_rate, double& rate_accumulator);
	int findEventPartner(const Event& next_event, const int& event_performer, const double& next_event_rate, double& rate_accumulator);
	void doFormationEvent(const int& event_performer, const int& event_partner);
	void doTradeEvent(const int& event_performer, const int& event_partner);
	void doRecEvent(const int& event_performer);
	void doNextEvent(const Event& next_event, const double& next_event_rate);
	
	void updateSupplyDemand();
	void removeSellers();
	void updateRiskScores();
	void resetDistinctSellers();
	
	void introduceDisease(int& init_inf);
	void introduceControl();
	
	void recalculateRates();
	double calcFormationRate(const int& buyer, const int& seller);
	double calcCessationRate(const int& buyer, const int& seller);
	double calcTradeRate(const int& buyer, const int& seller);
	double calcRecRate();
	
	double calcInfProb(const double& batch_size);
	double calcInfAnimals(const double& batch_size);
	bool canTransmit(const int& infectee, const int& infector);
	double calcDetectedAnimals(const double& inf_animals);
	
	double calcDemandRate(const int& farm);
	double calcSupplyRate(const int& farm);
	double calcPrice();
	
	bool isSeller(const int& buyer, const int& seller);
	bool isDistinctSeller(const int& buyer, const int& seller);
	
	void updateVector(std::vector<double>& vec, const double& update_amount);
	
	void runSim();
};

#endif