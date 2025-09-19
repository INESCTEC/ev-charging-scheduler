#include <iostream>
#include <fstream>  // For file handling
#include <vector>
#include <string>
#include <set>
#include "json.hpp"
#include "Highs.h"

// For convenience
using json = nlohmann::json;
using namespace std;


// charging forecasts for an EVSE
struct ChargingForecast {
    string evse_id;
    vector<double> power; // charging power values for the EVSE
};

// session data forecasts 
struct Session {
    string evse_id;
    string start_time;
    string end_time;
    double energy; // energy needed
    double capacity; // battery capacity
    vector<int> charging_index; // index of the timestamps when the vehicle is connected
    int charging_time; // number of timestamps when the vehicle is connected
    int evse; // EVSE index
};

//solar production forecasts
struct SolarForecasts{
    vector<string> timestamps;
    vector<double> values;
};



int main() {

    cout << endl << "Charging optimization" << endl;
   
    double freq = 12;

    
    double cost_hour = 10;
    double cost_energy = 10;
    double cost_soc = 2;
    double cost_dif = 0.001;

    vector<string> timestamps; // timestamps for the optimization window
    int total_charging_time = 0; // total connected time of all vehicles
    vector<string> unique_evse_ids; // EVSE IDs with charging sessions
    vector<ChargingForecast> chargingForecasts; // charging forecasts for each EVSE
    vector<vector<int>> stations; // EVSE indexes for each station
    SolarForecasts solar; // solar production forecasts
 
    // Read charging profiles forecasts
    string filePath = "/home/tulio/centralized_optimization/evse_forecasts/files/charging_forecasts.json";

    // Open the file
    ifstream inputFile(filePath);

    // Check if the file was opened successfully
    if (!inputFile.is_open()) {
        cerr << "Error: Could not open file " << filePath << endl;
        return 1; // Return an error code
    }

    // Parse the JSON data from the file
    json chargingData;
    inputFile >> chargingData;

    // Close the file
    inputFile.close();

    // Get timestamps for the optimization window
    const auto& firstEVSEArray = chargingData[0]; // First EVSE's array

    
    for (const auto& evseItem : firstEVSEArray["timestamps"]) {
        timestamps.push_back(evseItem["timestamp"]); // Store the timestamps
    }
    
    // Print the timestamps
    cout << "Timestamps: ";
    for (const auto& ts : timestamps) {
        cout << ts << " ";
    }
    cout << endl;



    // Read Session data forecasts
    std::string arrival_filename="/home/tulio/centralized_optimization/evse_forecasts/files/arrival_forecasts.json";
    std::ifstream input_file(arrival_filename);
    json arrivalData;
    input_file >> arrivalData;

    vector<Session> sessions;

    // For each session
    for (const auto& entry : arrivalData) {
        unique_evse_ids.push_back(entry["evse_id"]);
        for (const auto& session_data : entry["sessions"]) {
                    
            if (session_data["energy"] != 0 && ((timestamps.front() <= session_data["end_time"]) || (timestamps.back() <= session_data["start_time"]))){
                // if there is charging during the optimization window
                // Create a Session object and populate it with data
                Session session;
                session.evse_id = entry["evse_id"];
                session.start_time = session_data["start_time"];
                session.end_time = session_data["end_time"];
                session.energy = session_data["energy"];
                session.capacity = session_data["capacity"];

                sessions.push_back(session);
            }  
        }    
    }

    if (sessions.empty()){
        std::cerr << "Error: No sessions found" << std::endl;
        return 2;
    }


    // Count total charging time and charging indexes for each session
    for ( auto& session : sessions) {
        int charging_time = 0; //connected time of the vehicle
        vector<int> charging_index; // index of the timestamps when the vehicle is connected
        for (auto i = 0; i < timestamps.size(); ++i){
            if (timestamps[i] >= session.start_time && timestamps[i] <= session.end_time){
                // if the vehicle is connected at this timestamp
                charging_index.push_back(i);
                total_charging_time++;
                charging_time++;
            }
        }
        session.charging_index= charging_index;
        session.charging_time= charging_time;

    }


    int evse_count = unique_evse_ids.size();

    // Print the EVSE IDs
    cout << "Unique EVSE IDs:" << endl;
    for (const auto& evse_id : unique_evse_ids) {
        cout << evse_id << endl;
    }
    cout << endl << "-------------------" << endl;


    // Assign EVSE index to each session
    for (auto& session : sessions) {
        for (auto i = 0; i < evse_count; ++i) {
            if (unique_evse_ids[i]== session.evse_id)
                session.evse=i;
        }  
    } 

    // Print session data forecasts
    for (const auto& session : sessions) {
        cout << "EVSE ID: " << session.evse_id << endl;
        cout << "  Start: " << session.start_time << endl;
        cout << "  End: " << session.end_time << endl;
        cout << "  Energy: " << session.energy << endl;
        cout << "  Capacity: " << session.capacity << endl;
        cout << "  Charging time: " << session.charging_time << endl;
        cout << "  Evse number: " << session.evse << endl;
        for (const auto& index : session.charging_index) {
            cout << index << " ";
        }
        cout << endl << "-------------------" << endl;
    }

    

    // For each EVSE with a charging session
    for (auto i = 0; i < evse_count; ++i) {

        // Find charging data for the EVSE
        for (const auto& evseArray : chargingData) {
            if (unique_evse_ids[i]== evseArray["evse_id"]){

                ChargingForecast data;
                
                data.evse_id = evseArray["evse_id"]; 
                
                for (const auto& evseItem : evseArray["timestamps"]) 
                    data.power.push_back(evseItem["power"]);      // Save power value

                // Add the populated ChargingForecast to the list
                chargingForecasts.push_back(data);
            }           
        }
    }

  

    // Print charging forecasts for each EVSE
    for (const auto& data : chargingForecasts) {
        cout << "EVSE ID: " << data.evse_id << endl;
        cout << "Power values: ";
        for (const auto& pwr : data.power) {
            cout << pwr << " ";
        }
        cout << endl;
    }
    



    // Read charging stations data
    std::ifstream stations_file("/home/tulio/centralized_optimization/evse_forecasts/files/stations.json");

    json station_data;
    stations_file >> station_data;

    

    // For each station
    for (const auto& entry : station_data) {
        vector<int> evses;
        for( auto &array : entry["evse_ids"]){
            string evse_id=array["evse_id"].get<std::string>();
            for (auto i = 0; i < evse_count; ++i) {
                if (unique_evse_ids[i]== evse_id)
                    evses.push_back(i); // Save the EVSE index
            }  
        } 
        if (!evses.empty()) {
            stations.push_back(evses); // Add the EVSE indexes to the stations vector
        }
        
    }

    // // Print the stations and their EVSEs
    // for (const auto& station : stations) {
    //     for (const auto& evse : station) {
    //         cout << evse << " ";
    //     }
    //     cout << endl << "-------------------" << endl;
    // }


    // Read solar production forecasts
    std::ifstream solar_file("/home/tulio/centralized_optimization/evse_forecasts/files/solar_forecast.json");
    json solar_data;
    solar_file >> solar_data;


    for (const auto& forecast : solar_data.at("forecast")) {
        // Find first timestamp in the optimization window
        if (forecast.at("timestamp").get<std::string>() >= timestamps[0]){
            for (int i = 0; i < 3; ++i) {
                solar.timestamps.push_back(forecast.at("timestamp").get<std::string>());
                solar.values.push_back(forecast.at("value").get<double>());
            }
        }
        
    }

    // for (size_t i = 0; i < solar.timestamps.size(); ++i) {
    //     std::cout << "Timestamp: " << solar.timestamps[i] << ", Value: " << solar.values[i] << std::endl;
    // }


    



    HighsModel model;
    model.lp_.sense_ = ObjSense::kMinimize; // Set the objective to minimize

    // Number of variables
    int power_var = timestamps.size()*evse_count; // charging power for each EVSE at each timestamp
    int soc_var = sessions.size(); // absolute deviation of the SOC for each EVSE
    int charging_diff_var = total_charging_time; // absolute deviation of the charging power for each EVSE
    int variable_n= power_var + soc_var + charging_diff_var; // total number of variables


    // Number of constraints
    int station_const = timestamps.size()*stations.size(); // power limits for each charging station at each timestamp
    int energy_const = sessions.size(); // charging limited to energy need
    int soc_const = 2* sessions.size(); // absolute deviation of the SOC for each EVSE
    int charging_diff_const = 2 * total_charging_time; // absolute deviation of the charging power for each EVSE
    int constraint_n= station_const + energy_const + soc_const + charging_diff_const; // total number of constraints


    cout << "Variable number: " <<variable_n<< endl;
    cout << "Constraint number: " <<constraint_n<< endl;


    model.lp_.num_col_ = variable_n;
    model.lp_.num_row_ = constraint_n;


    // Setting variable bounds
	std::vector<double> lb(variable_n, 0);
	std::vector<double> ub(variable_n, 0);


    // Charging power limited to EVSE maximum power, when there is a session, 0 otherwise
    for (auto i = 0; i < sessions.size(); ++i) {
        for (const auto&j : sessions[i].charging_index){
            ub[j + sessions[i].evse * timestamps.size()] = 7360;
        }
    }

    
    // other variables are not limited
    for (auto i = power_var; i < variable_n; ++i) {
		lb[i]=-1.0e30;
        ub[i]=1.0e30;
	}

	model.lp_.col_lower_ = lb;
	model.lp_.col_upper_ = ub;

    // // Print variable bounds
    // cout << "Variable bounds" << endl;
    // for (const auto& l : lb) {
    //     cout << l << " ";
    // }
    // cout << endl ;
    // for (const auto& l : ub) {
    //     cout << l << " ";
    // }
    // cout << endl << "-------------------" << endl;




    // Setting constraint bounds
	std::vector<double> L(constraint_n, -1.0e30);
	std::vector<double> U(constraint_n, 0);

    double average_soc_need=0; // average of the SOC needed for all vehicles
    for (auto i = 0; i < sessions.size(); ++i) {
		average_soc_need += static_cast<double>(sessions[i].energy)/(sessions[i].capacity*sessions.size());

	}

    // constraint for the absolute deviation of the SOC
	for (auto i = 0; i < sessions.size(); ++i) {
		L[i + station_const + energy_const + soc_const/2] = (static_cast<double>(sessions[i].energy)/sessions[i].capacity) - average_soc_need; // SOC deviation for that EVSE
	}

    // constraint for the absolute deviation of the charging power
    int x=0;
    for (auto i = 0; i < sessions.size(); ++i) {
        for (auto j = 0; j < sessions[i].charging_time; ++j) {
            L[x+constraint_n-charging_diff_const/2]=chargingForecasts[sessions[i].evse].power[sessions[i].charging_index[j]]; //forecasted charging power for that EVSE/timestamp
            x++;
        }
    }

    // Power limits for each charging station at each timestamp
    for (auto i = 0; i < timestamps.size()*stations.size(); ++i) {
		U[i] = 10000;
	}

    // Charging limited to energy need
    for (auto i = 0; i < sessions.size(); ++i) {
		U[station_const+i] = sessions[i].energy * freq;
	}

    // constraint for the absolute deviation of the SOC
    for (auto i = 0; i < sessions.size(); ++i) {
		U[i+station_const+ energy_const] = (static_cast<double>(sessions[i].energy)/sessions[i].capacity) - average_soc_need;  // SOC deviation for that EVSE
	}
    
    // Constraints with no upper limit
    for (auto i = station_const + energy_const + soc_const/2; i < constraint_n; ++i) {
		U[i] = 1.0e30;
	}

    // constraint for the absolute deviation of the charging power
    x=0;
    for (auto i = 0; i < sessions.size(); ++i) {
        for (auto j = 0; j < sessions[i].charging_time; ++j) {
            U[x+constraint_n- charging_diff_const]=chargingForecasts[sessions[i].evse].power[sessions[i].charging_index[j]]; //forecasted charging power for that EVSE/timestamp
            x++;
        }
    }

    // // Print constraint bounds
    // cout << "Constraint bounds" << endl;
    // for (const auto& l : L) {
    //     cout << l << " ";
    // }
    // cout << endl ;
    // for (const auto& l : U) {
    //     cout << l << " ";
    // }
    // cout << endl << "-------------------" << endl;


	model.lp_.row_lower_ = L;
	model.lp_.row_upper_ = U;





    // Constructing the constraint matrix

    // Matrix size:  
    int a_station= timestamps.size()*evse_count; // charging power for each EVSE at each timestamp
    int a_energy= total_charging_time; // charging powers for each EVSE when there is a session
    int a_soc= 2* sessions.size()* (total_charging_time+1); // for each EVSE, absolute deviation of the SOC and charging powers of all EVSE when there is a session
    int a_charging_diff= 4 * total_charging_time; // or each EVSE when there is a session, charging power and its absolute deviation

    int a_size= a_station + a_energy + a_soc + a_charging_diff; // total number of variables in the constraint matrix


    std::vector<double> a(a_size, 1); // constraint matrix
	std::vector<int> index(a_size, 0); // index of the variables in the constraint matrix
	std::vector<int> start(constraint_n + 1, 0); // start of each constraint in the matrix

    cout << "A size: " <<a_size<< endl;

    // constraints for the absolute deviation of the SOC, cost of the charging power for all sessions
    x=a_station+a_energy;
    for (auto i = 0; i < sessions.size(); ++i) {
        for (auto j = 0; j < sessions.size(); ++j) {
            for (auto k = 0; k < sessions[j].charging_time; ++k) {
                a[x]=-1/(freq*sessions[j].energy*sessions.size());
                a[x+a_soc/2]=-1/(freq*sessions[j].energy*sessions.size());
                x++;

            }
        }
        x++;
	}

    // constraints for the absolute deviation of the SOC, cost of the charging power for the session 
    x=a_station+a_energy;
    for (auto i = 0; i < sessions.size(); ++i) {
        for (auto j = 0; j < sessions[i].charging_time; ++j) {
            a[x]+=1/(freq*sessions[i].energy);
            a[x+a_soc/2]+=1/(freq*sessions[i].energy);
            x++;
        }
        x += total_charging_time+1;
    }

    // absolute deviation of the SOC upper limit
    for (auto i = 0; i < sessions.size(); ++i) {
        a[a_station+a_energy+total_charging_time+i*(total_charging_time+1)]=-1;
    }

    // absolute deviation of the charging power upper limit
    for (auto i = 0; i < total_charging_time; ++i) {
        a[a_size-a_charging_diff+1 +2*i]=-1;
    }

    // // Print the A matrix
    // cout << "A: " << endl;
    // for (const auto& l : a) {
    //     cout << l << " ";
    // }


    // constraints for charging station power limits
    x=0;
    for (auto i = 0; i < timestamps.size(); ++i) {
        for (const auto& station : stations) {
            for (const auto& evse_id : station) {
                index[x] = i+timestamps.size()*evse_id;
                x++;
            }
        }
	}
    // cout << endl<< x;

    // constraints for charging limited to energy need
    for (auto i = 0; i < sessions.size(); ++i) {
        for (const auto&j : sessions[i].charging_index){
            index[x] = j + sessions[i].evse * timestamps.size();
            x++;
        }
    }
    // cout << endl<< x <<endl;

    // index of the charging powers in the constraint for the absolute deviation of the SOC
    for (auto i = 0; i < sessions.size(); ++i) {
        for (const auto&j : sessions[i].charging_index){
            for (auto k = 0; k < 2*sessions.size(); k++){
                index[x+k*(total_charging_time+1)] = j + sessions[i].evse * timestamps.size();       
            }
            x++;
        }
    }
    
    // index of the absolute deviation of the SOC
        x=a_station+a_energy+total_charging_time;
    for (auto k = 0; k < sessions.size(); k++){
        index[x+k*(total_charging_time+1)] =  power_var+k; // upper limit constraint
        index[x+(k+sessions.size())*(total_charging_time+1)] =  power_var+k; // lower limit constraint
        
    }

    // index of the charging powers in the constraint for the absolute deviation of the charging power
    x=a_size-a_charging_diff;
    for (auto i = 0; i < sessions.size(); ++i) {
        for (const auto&j : sessions[i].charging_index){
            index[x] = j + sessions[i].evse * timestamps.size(); // upper limit constraint
            index[x+a_charging_diff/2] = j + sessions[i].evse * timestamps.size(); // lower limit constraint
            x+=2;
        }
    }

    // index of the absolute deviation of the charging power
    x=a_size-a_charging_diff+1;
    int var=power_var+soc_var;
    for (auto i = 0; i < sessions.size(); ++i) {
        for (const auto&j : sessions[i].charging_index){
            index[x] = var; // upper limit constraint
            index[x+a_charging_diff/2] = var; // lower limit constraint
            x+=2;
            var++;
        }
    }

    // // Print the index matrix
    // cout << endl << "Index: " << endl;
    // for (const auto& l : index) {
    //     cout << l << " ";
    // }
    


    // constraints for charging station power limits
    for (auto x = 0; x < timestamps.size(); ++x) {
        for (auto i = 0; i < stations.size(); ++i) {
            start[i+1+x*stations.size()]=start[i+x*stations.size()]+stations[i].size(); //numbrer of EVSEs in the station
        }
    }

    // constraints for charging limited to energy need
    x=station_const+1;
    for (auto i = 0 ; i <  energy_const; ++i) {
        start[i+x]=start[i+x-1]+sessions[i].charging_time; // number of timestamps when the vehicle is connected

    }
    
    // constraints for the absolute deviation of the SOC
    for (auto i = station_const+energy_const+1 ; i < station_const+energy_const+1+soc_const; ++i) {
        start[i]=start[i-1]+(total_charging_time+1); // number of timestamps when vehicles are connected + absolute deviation of the SOC
    }

    // constraints for the absolute deviation of the charging power
    for (auto i = station_const+energy_const+soc_const+1 ; i < constraint_n+1; ++i) {
        start[i]=start[i-1]+2; // charging power + absolute deviation of the charging power
    }

    // // Print the start matrix
    // cout << endl << "Start: " << endl;
    // for (const auto& l : start) {
    //     cout << l << " ";
    // }


    model.lp_.a_matrix_.start_ = start;
	model.lp_.a_matrix_.index_ = index;
	model.lp_.a_matrix_.value_ = a;

    model.lp_.a_matrix_.format_ = MatrixFormat::kRowwise;

    // Costs for the objective function
    std::vector<double> cost(variable_n, 1);

    cout << endl << "Solar: " << endl;
    
    // Cost coefficient for transfered energy
    for (auto i = 0 ; i < power_var; ++i) {
        cost[i]=-cost_energy;
    }

    // Cost coefficient for the solar production (added to the same variables)
    for (auto i = 0; i < evse_count; ++i) {
        for (auto j = 0; j < timestamps.size(); ++j) {
            cost[j+i*timestamps.size()]+=cost_hour*(500/(solar.values[j]+1));
            cout << solar.values[j] << " ";
        }
    }
    cout << endl;
    // Cost coefficient for the absolute deviation of the SOC
    for (auto i = power_var ; i < power_var+soc_var; ++i) {
        cost[i]=cost_soc;
    }
    
    // Cost coefficient for the absolute deviation of the charging power
    for (auto i = power_var+soc_var ; i < variable_n; ++i) {
        cost[i]=cost_dif;
    }
    
    // // Print the cost matrix
    // cout << endl << "Cost: " << endl;
    // for (const auto& l : cost) {
    //     cout << l << " ";
    // }
    // cout << endl;

    model.lp_.col_cost_ = cost;

    // Create a Highs instance
    Highs highs;
    HighsStatus return_status;
    
    // Pass the model to HiGHS
    return_status = highs.passModel(model);
  
    // Get a const reference to the LP data in HiGHS
    const HighsLp& lp = highs.getLp();
    
    // Solve the model
    return_status = highs.run();
    assert(return_status==HighsStatus::kOk);
    
    // Get the model status
    const HighsModelStatus& model_status = highs.getModelStatus();
    // Check if the model status is optimal
    if (model_status != HighsModelStatus::kOptimal) {
        std::cerr << "Error: Model did not reach optimal status." << std::endl;
        return 0;
    }

    // Print model information
    const HighsInfo& info = highs.getInfo();
    cout << "Simplex iteration count: " << info.simplex_iteration_count << endl;
    cout << "Objective function value: " << info.objective_function_value << endl;
    cout << "Primal  solution status: " << highs.solutionStatusToString(info.primal_solution_status) << endl;
    cout << "Dual    solution status: " << highs.solutionStatusToString(info.dual_solution_status) << endl;
    cout << "Basis: " << highs.basisValidityToString(info.basis_validity) << endl;


    const bool has_values = info.primal_solution_status;
    const HighsSolution& solution = highs.getSolution();

    // If there is a solution
    if (has_values) {

        // Print the solution
        cout << endl << "Results: " << endl;
		for (auto i = 0; i < evse_count; ++i)
		{
            for (auto j = 0; j < timestamps.size(); ++j){
               cout << solution.col_value[i*timestamps.size()+j] << " "; 
            }
			cout << endl;
		}
        cout << endl;
        // for (auto i = power_var; i < power_var+soc_var; ++i)
		// {
		// 	cout << solution.col_value[i] << " ";
		// }
        // cout << endl;
        // for (auto i = power_var+soc_var; i < variable_n; ++i)
		// {
		// 	cout << solution.col_value[i] << " ";
		// }


        // Save the results to a JSON file
        json json_content;

        std::string outputFileName = "/home/tulio/centralized_optimization/optimization/files/optimization_results.json";

        std::ofstream o( outputFileName );
        json evses_json = json::array();

        // For each EVSE
        for (auto i = 0; i < evse_count; ++i)
        {
            json new_evse;

            new_evse["evse_id"] = chargingForecasts[i].evse_id;
            json measurements = json::array();
            for (auto j = 0; j < timestamps.size(); ++j){
                
                json new_measurement;
                new_measurement["timestamp"] = timestamps[j];
                new_measurement["power"] = solution.col_value[i*timestamps.size()+j];
                measurements.push_back( new_measurement );
            }
            new_evse["measurements"] = measurements;
            evses_json.push_back( new_evse );
        }
        json_content["optimization"] = evses_json;

        o << std::setw( 4 ) << json_content << std::endl;

        o.close();
        cout << endl << outputFileName;

    
	}


    return 0;
}