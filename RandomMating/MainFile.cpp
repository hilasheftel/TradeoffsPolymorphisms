#include <algorithm> 
#include <assert.h>
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>
#include <chrono>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <numeric>
#include "person.h"
#include <stdio.h>
#include <time.h>
#include <math.h>       
#include <list>
#include <regex>

using namespace std;
using namespace std::placeholders;
using namespace std::chrono;
default_random_engine generator;
common_person common_pers;
double PI = 3.14159265;
double archetype;
double theta1;
double theta2;
double lambda;

const std::string currentDateTime() {
	time_t     now = time(0);
	struct tm  tstruct;
	char       buf[80];
	tstruct = *localtime(&now);

	strftime(buf, sizeof(buf), "_%Y_%m_%d_%H_%M_%S", &tstruct);

	return buf;
}

milliseconds currentDateTimeMili() {
	milliseconds ms = duration_cast< milliseconds >(
		system_clock::now().time_since_epoch()
		);
	return ms;
}


nanoseconds currentDateTimeNano() {
	nanoseconds ns = duration_cast< nanoseconds >(
		system_clock::now().time_since_epoch()
		);
	return ns;
}

double my_calc_perf(double x, double y, double x0, double theta){
	double ans = ((((x - x0)*(x - x0)) + (y*y*lambda*lambda))*cos(theta)*cos(theta)) + (((y*y) + ((x - x0)*(x - x0)*lambda*lambda))*sin(theta)*sin(theta)) - ((x - x0)*y*(-1 + (lambda*lambda))*sin(2 * theta));
	return ans;
}

void writePareto(vector<person>& pop, string filename, int gen) {
	ofstream file2Write;
	ostringstream filename2;
	filename2 << filename << "_currentgen_" << gen << ".tsv";
	//filename2 << filename << "_currentgen_" << gen << ".tsv";
	//filename2 << "D:\\debug\\" << filename << "_currentgen_" << gen << ".tsv";

	//filename2 << "D:\\SimulationResults\\" << filename << "_currentgen_" << gen << ".tsv";
	//filename2 << "C:\\Users\\Alonlab\\Copy\\Simulations\\" << filename;
	//results_SimpleFitness_popsize_" << popSize << "_tournament_" << tournamentSize <<"_rec_" <<recombRate <<"_mut_" <<  mutRate<<"_beta_" << selStrength<<".tsv";
	file2Write.open(filename2.str());
	for (unsigned int i = 0; i < pop.size(); i++) {
		string a = pop.at(i).print();
		file2Write << a;
	}
	file2Write.close();

}

void writeParetofixed(string filename, int gen) {
	ofstream file2Write;
	ostringstream filename2;
	filename2 << filename << "_fixed_currentgen_" << gen << ".tsv";
	//filename2 << filename << "_currentgen_" << gen << ".tsv";
	//filename2 << "D:\\debug\\" << filename << "_currentgen_" << gen << ".tsv";

	//filename2 << "D:\\SimulationResults\\" << filename << "_currentgen_" << gen << ".tsv";
	//filename2 << "C:\\Users\\Alonlab\\Copy\\Simulations\\" << filename;
	//results_SimpleFitness_popsize_" << popSize << "_tournament_" << tournamentSize <<"_rec_" <<recombRate <<"_mut_" <<  mutRate<<"_beta_" << selStrength<<".tsv";
	file2Write.open(filename2.str());
	string a = common_pers.print();
	file2Write << a;

	file2Write.close();

}



bool compare_fitnesses(double f1, double f2, vector<double>& data)
{

	return data.at(f1) < data.at(f2);
}


bool compare_p1s(double f1, double f2, vector<double>& p1s, vector<double>& p2s)
{
	double numerror = 1e-30;
	if (abs(p1s.at(f1) - p1s.at(f2))< numerror){
		return p2s.at(f1) < p2s.at(f2);
	}
	return p1s.at(f1) < p1s.at(f2);
}
bool compare_p1s_pop(person& p1, person& p2){
	return p1.perf1 < p2.perf1;
}





void shuffle_element(vector<int>& vec, int ind1, int ind2){

	int temp = vec.at(ind1);
	vec.at(ind1) = vec.at(ind2);
	vec.at(ind2) = temp;

}

bool compare_first(pair<double, int> i, pair<double, int> j) {
	return i.first <j.first;
}
vector <person> Biological2(vector <person>& pop, int popSize, int gen){
	//mark every junction with names and values. if cut above check and add according to multiplicity if cuts below same

	vector <int> alphas(popSize, 0);
	std::iota(std::begin(alphas), std::end(alphas), 0);

	//vector <bool> persons(popSize, false); //is a person there to begin the layer

	vector<person> offsprings;
	std::shuffle(alphas.begin(), alphas.end(), generator);

	for (unsigned int i = 0; i < popSize; i++)
	{

		int selected_alpha = alphas.at(i);

		vector <pair<double, int>> fitness(pop.size());

		for (unsigned int j = 0; j < pop.size(); j++){
			double alpha_normalized = (float(selected_alpha) / (popSize - 1));

			fitness.at(j) = pair<double, int>((1 - alpha_normalized)* pop.at(j).perf1 + alpha_normalized* pop.at(j).perf2, j);
		}
		pair<double, int> smallest = *std::min_element(fitness.begin(), fitness.end(), compare_first);
		int best_at_niche = smallest.second;
		person selected_person(pop.at(best_at_niche));
		selected_person.set_niche(selected_alpha);
		offsprings.push_back(person(selected_person));
		pop.erase(pop.begin() + best_at_niche);

	}





	return offsprings;

}


string find_parameter(string filename, string parameter_name){
	regex re(parameter_name);

	std::smatch m;
	std::smatch m2;

	std::regex_search(filename, m, re);
	regex dig("[0-9]+(_|.)");
	string to_match = string(m.suffix());
	std::regex_search(to_match, m2, dig);
	string match = m2.str();
	match.erase(match.end() - 1);
	return match;
}


void simulation(int GENS, int popSize, double mutRate, string filename, double initial_mean, double initial_x, double initial_y, int selection_strength, bool test, int save_every, double sigma, double top, string input_file) {
	//unordered_map <double, int> fixed_mutations;
	int numberOfMutations;
	uniform_int_distribution<int> randfilename(0, 1000);
	filename += currentDateTime();
	filename += "_rand_";
	filename += std::to_string(randfilename(generator));
	vector <person> pop;
	int start_gen = -1;

	uniform_real_distribution<double> dist(0.0, 1.0);
	double initial_mutation_r = dist(generator);
	if (input_file == "no_input_file"){
		for (int i = 0; i < popSize; i++) {

			//	person per();
			pop.push_back(person(i));
			pop.at(i).add_polymorphism(initial_mutation_r, initial_x, initial_y);
			pop.at(i).mutate(initial_mean);
		}
	}
	else{
		string input_fixed_file = string(input_file);
		std::size_t found = input_fixed_file.find("_currentgen");
		input_fixed_file.insert(found, "_fixed");

		//regex re_prefix("_currentgen");
		//std::smatch m_prefix;
		//bool res = std::regex_search(input_file, m_prefix, re_prefix);
		//filename = string(m_prefix.prefix());
		filename = input_file.substr(0, found);
		ifstream infile;
		//		int pop_size = stoi(find_parameter(input_file, "pop"));
		//		int currentgen = ;
		std::size_t found_point = input_file.find_last_of(".");
		string current_gen = input_file.substr(found + 12, found_point - (found + 12));
		start_gen += stoi(current_gen);

		common_pers.reset_common_person(input_fixed_file);
		//common_pers.set_pop_size(pop_size);
		infile.open(input_file);
		string line;
		int kk = 0;
		while (!infile.eof()){
			std::getline(infile, line);
			if (line == ""){
				break;
			}
			kk++;
			pop.push_back(person(line));

		}

	}
	//Evolution begins! 



	for (int gen = start_gen + 1; gen < GENS; gen++) {

		if (gen % save_every == 0){

			//if (gen==10 || gen==100000) {
			writePareto(pop, filename, gen);
			writeParetofixed(filename, gen);
		}
		vector <person> offsprings;
		common_pers.reset();


		for (int i = 0; i< selection_strength; i++) {
			uniform_int_distribution<int> choose_parent(0, pop.size() - 1);

			// recombination			
			//this code should be used in case recombination does not take into accoung linkage between genes
			int parent1 = choose_parent(generator);
			int parent2 = choose_parent(generator);
			while (parent2 == parent1){
				parent2 = choose_parent(generator);
			}
			Chromosome g1 = pop.at(parent1).ForInfiniteRecombination(mutRate);
			Chromosome g2 = pop.at(parent2).ForInfiniteRecombination(mutRate);

			person p(i, g1, g2);
			offsprings.push_back(p);
		}

		//Start no tour
		pop = offsprings;
		/*if ((gen + 1) % save_every == 0){

		//if (gen==10 || gen==100000) {
		writePareto(pop, filename + "_chosen", gen + 1);
		writeParetofixed(filename + "_chosen", gen + 1);
		}*/
		//milliseconds ms1 = currentDateTimeMili();
		pop = Biological2(pop, popSize, gen);
		//milliseconds ms2 = currentDateTimeMili();
		/*if (gen % 100 == 0){
		std::cout << "biological took at gen : "<< gen
		<<"\n" << (ms2 - ms1).count()
		<< " milseconds\n";
		}*/

		//<< " miliseconds\n";
		// kids are thought upon.
		// newPareto(pop);
		//		milliseconds ms1 = currentDateTimeMili();

		//	milliseconds ms2 = currentDateTimeMili();




	}
	//The end of the world - now collect data 


}
int main(int argc, char* argv[])  {
	//get 9 parameters - simple fitness: Num generations, pop size, mut rate, file name, alpha, initial mutation mean,  lambda, archetype, repeats

	if (argc < 13) {
		std::cerr << "Usage: " << argv[0] << "SOURCE DESTINATION" << std::endl;
		cout << "the right order of the parameters is:\n1. number of generations\n, 2. population size\n 3. mutation rate\n 4. file name template\n 5. alpha\n 6. initial mutation number mean\n 7. lambda\n 8. archetype\n 9. repeats\n 10. initial x\n 11. initial y\n 12. num layers\n 13. highest survival\n14.testmode\n 15. num bins \n16. binned selection \n17. initiate from file \n18. filename!";
		return 11;
	}

	std::cout.flush();

	chrono::high_resolution_clock myclock;
	chrono::high_resolution_clock::time_point beg = myclock.from_time_t(rand());
	chrono::high_resolution_clock::duration d = myclock.now() - beg;

	unsigned seed2 = d.count();
	generator.seed(seed2);
	// the parameters: 1. number of generations, 2. population size 3. mutation rate 4. file name template 6. initial mutation number mean 7. lambda 8. archetype 9. repeats 10. initial x 11. initial y
	const int GENS = atoi(argv[1]);
	const int popSize = atoi(argv[2]); //tournamentSize = atoi(argv[2]);
	double mutRate = atof(argv[3]);
	string filename = argv[4];
	double initial_mean = atof(argv[5]);
	lambda = atof(argv[6]);
	archetype = atof(argv[7]);
	int repeats = atoi(argv[8]);
	double initial_x = atof(argv[9]);
	double initial_y = atof(argv[10]);
	double selection_factor = atof(argv[11]);
	double sigma = atof(argv[12]);
	bool test = bool(atoi(argv[13]));
	int save_every = atoi(argv[14]);
	theta1 = atof(argv[15]);
	theta2 = atof(argv[16]);
	double top = atof(argv[17]);



	std::string final_res;
	string temp;
	int dotIndex;
	std::cout << final_res << std::endl;
	filename += "_pop_";
	filename += to_string(popSize);
	filename += "_mut_";
	temp = to_string(mutRate);
	dotIndex = temp.find(".");
	final_res = temp.substr(0, dotIndex) + "dot" + temp.substr(dotIndex + 1, 3);
	filename += final_res;
	filename += "_initial_";
	temp = to_string(initial_mean);
	dotIndex = temp.find(".");
	final_res = temp.substr(0, dotIndex);
	filename += final_res;
	filename += "_lambda_";
	temp = to_string(lambda);
	dotIndex = temp.find(".");
	final_res = temp.substr(0, dotIndex) + "dot" + temp.substr(dotIndex + 1, 3);
	filename += final_res;
	filename += "_arc_";
	temp = to_string(archetype);
	dotIndex = temp.find(".");
	final_res = temp.substr(0, dotIndex);
	filename += final_res;

	filename += "_select_";
	temp = to_string(selection_factor);
	dotIndex = temp.find(".");
	final_res = temp.substr(0, dotIndex) + "dot" + temp.substr(dotIndex + 1, 1);
	filename += final_res;
	filename += "q1";

	temp = to_string(theta1);
	dotIndex = temp.find(".");
	final_res = temp.substr(0, dotIndex) + "dot" + temp.substr(dotIndex + 1, 2);
	filename += final_res;

	filename += "q2";

	temp = to_string(theta2);
	dotIndex = temp.find(".");
	final_res = temp.substr(0, dotIndex) + "dot" + temp.substr(dotIndex + 1, 2);
	filename += final_res;

	filename += "_top_";

	temp = to_string(top);
	dotIndex = temp.find(".");
	final_res = temp.substr(0, dotIndex) + "dot" + temp.substr(dotIndex + 1, 2);
	filename += final_res;

	int selection_strength = floor(selection_factor*popSize);

	bool from_file = bool(atoi(argv[18]));
	string input_file;
	if (from_file){
		input_file = argv[19];
	}
	else{
		input_file = "no_input_file";
	}
	common_pers.set_pop_size(selection_strength);
	if (test) {
		vector <person> pop;
		vector <person> pop_2;
		vector <person> pop_3;

		pop.push_back(person(1));
		pop.at(0).set_perf_1(6.789);
		pop.at(0).set_perf_2(5.223);
		pop.at(0).set_phen_x(1);
		pop.at(0).set_phen_y(1);
		pop.push_back(person(1));
		pop.at(1).set_perf_1(0.21);
		pop.at(1).set_perf_2(1.34);
		pop.at(1).set_phen_x(2);
		pop.at(1).set_phen_y(2);
		pop.push_back(person(1));
		pop.at(2).set_perf_1(0.82);
		pop.at(2).set_perf_2(1.11);
		pop.at(2).set_phen_x(3);
		pop.at(2).set_phen_y(3);
		pop.push_back(person(1));
		pop.at(3).set_perf_1(0.63);
		pop.at(3).set_perf_2(1.21);
		pop.at(3).set_phen_x(4);
		pop.at(3).set_phen_y(4);
		pop.push_back(person(1));
		pop.at(4).set_perf_1(4.32);
		pop.at(4).set_perf_2(6.37);
		pop.at(4).set_phen_x(5);
		pop.at(4).set_phen_y(5);
		pop.push_back(person(1));
		pop.at(5).set_perf_1(0.21);
		pop.at(5).set_perf_2(1.34);
		pop.at(5).set_phen_x(2);
		pop.at(5).set_phen_y(2);
		//		Biological(pop, 6, 1);
	}

	else {
		for (int j = 0; j < repeats; j++) {
			simulation(GENS, popSize, mutRate, filename, initial_mean, initial_x, initial_y, selection_strength, test, save_every, sigma, top, input_file);
		}

	}
	return 0;
}
