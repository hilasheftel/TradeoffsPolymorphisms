#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <random>
#include <algorithm>    // std::sort
#include "person.h"
#include <math.h>       


extern default_random_engine generator;
extern common_person common_pers;
extern double archetype;
extern double theta1;
extern double theta2;
extern double lambda;
extern double PI;

double calc_perf(double x, double y, double x0, double theta){
	double ans = ((((x - x0)*(x - x0)) + (y*y*lambda*lambda))*cos(theta)*cos(theta)) + (((y*y) + ((x - x0)*(x - x0)*lambda*lambda))*sin(theta)*sin(theta)) - ((x - x0)*y*(-1 + (lambda*lambda))*sin(2 * theta));
	return ans;
}

person::person(int bornniche) {
	phen_x = 0.0;
	phen_y = 0.0;
	perf1 = 0.0;
	perf2 = 0.0;
	niche = bornniche;
}

person::person(int bornniche, Chromosome a, Chromosome b) {
	gen1 = a;
	gen2 = b;
	pair<double, double> phen = phenotype();
	phen_x = phen.first;
	phen_y = phen.second;
	perf1 = calc_perf(phen_x, phen_y, 0.0, theta1);
	perf2 = calc_perf(phen_x, phen_y, archetype, theta2);
	niche = bornniche;

}

size_t hash_combine(size_t lhs, size_t rhs) {
	lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
	return lhs;
}

person::person(string str){
	//ifstream infile;
	//infile.open(filename);
	//string 
	//while (!infile.eof){
	std::vector<std::string> tokens;
	std::istringstream iss(str);
	std::string token;
	bool isgen1 = false;
	bool isgen2 = false;
	int order = 0;
	double r;
	double x;
	double y;
	while (std::getline(iss, token, '\t')) 
		tokens.push_back(token);
	
	for (int i = 0; i < tokens.size(); i++){
		if (tokens.at(i) == "gen1: "){
			isgen1 = true;
			continue;
		}
		if (tokens.at(i) == "gen2: "){
			isgen1 = false;
			isgen2 = true;
			continue;
		}
		switch (order){
			case 0: //r
				r = stod(tokens.at(i));
				order += 1;
				break;
			case 1://x
				x = stod(tokens.at(i));
				order += 1;
				break;
			case 2://y
				y = stod(tokens.at(i));
				order = 0;
				Gene temp(r, x, y);
				if (isgen1){
					gen1.push_back(temp);
				}
				else{
					gen2.push_back(temp);
				}
				break;
		}
	}
	pair<double, double> phen = phenotype();
	phen_x = phen.first;
	phen_y = phen.second;
	perf1 = calc_perf(phen_x, phen_y, 0.0, theta1);
	perf2 = calc_perf(phen_x, phen_y, archetype, theta2);

}
Chromosome person::get_gene(int geneNum)
{
    if (geneNum == 0)
        return gen1;
    else
        return gen2;

}


bool sort_genome(const Gene& t1, const Gene& t2){
    return t1.r < t2.r;
}

Chromosome person::ForRecombination(double mutRate) {
	uniform_real_distribution<double> recomb_split(0.0, 1.0);

	double temp = recomb_split(generator);
	Chromosome res;
	for (unsigned int i = 0; i < gen1.size(); i++) {
		if (gen1.at(i).r < temp) {
			res.push_back(gen1.at(i));
		}
	}

	for (unsigned int i = 0; i < gen2.size(); i++) {
		if (gen2.at(i).r >= temp) {
			res.push_back(gen2.at(i));
		}
	}
	uniform_real_distribution<double> dist(0.0, 1.0);
	uniform_real_distribution<double> angle(0.0, 6.283185);
	//std::exponential_distribution<double> effectSize(1);
	std::gamma_distribution<double> effectSize(1, 1);
	//uniform_real_distribution<double> epistasis(-epistasis_var, epistasis_var);
	uniform_int_distribution<int> DistZeroOne(0, 1);
	poisson_distribution<int> mutNum(mutRate);

	int mutation_num = mutNum(generator);
	for (int i = 0; i < mutation_num; i++) {

		double mut = dist(generator);
		double ang = angle(generator);
		//double eff = 1;
		double eff = sqrt(effectSize(generator));
		//double ep = epistasis(generator);
		double ep_angle = angle(generator);

		//double epx = ep*cos(ap_angle);
		//double epy = ep*sin(ap_angle);

		//Gene temp(mut, eff*cos(ang), eff*sin(ang), epx, epy);
		Gene temp(mut, eff*cos(ang), eff*sin(ang));
		if (!(common_pers.increase_gene(temp))){
			res.push_back(temp);
		}
	}
	return res;
}

Chromosome person::ForInfiniteRecombination(double mutRate) {
    uniform_int_distribution<int> intDistZeroOne (0,1);
    Chromosome res;
    Chromosome comb(gen1);
    comb.insert(comb.end(), gen2.begin(), gen2.end());
    sort(comb.begin(), comb.end(), sort_genome);

    for (unsigned int i = 0 ; i < comb.size(); i++) {
		Gene g = comb.at(i);
		if (common_pers.in_fixed(g)){
			continue;
		}
		
        if (i < comb.size()-1 && comb.at(i)==comb.at(i+1)){
			if (!(common_pers.increase_gene(g))){
				res.push_back(g);
			}
            i++;
            continue;
        }

        if(intDistZeroOne(generator)==1) {
			if (!(common_pers.increase_gene(g))){
				res.push_back(g);
			}

        }
    }

	uniform_real_distribution<double> dist(0.0,1.0);
    uniform_real_distribution<double> angle(0.0,6.283185);
    //std::exponential_distribution<double> effectSize(1);
	std::gamma_distribution<double> effectSize(1, 1);
	uniform_int_distribution<int> DistZeroOne (0,1);
	poisson_distribution<int> mutNum(mutRate);

	int mutation_num = mutNum(generator);
	for (int i = 0; i < mutation_num; i++) {
 
		double mut =  dist(generator);
		double ang = angle(generator);
		//double eff = 1;
		double eff = sqrt(effectSize(generator));
		//double choose = DistZeroOne(generator);
		//double eff = eff1*choose+eff2*(1-choose);

		Gene temp(mut, eff*cos(ang), eff*sin(ang));
		if (!(common_pers.increase_gene(temp))){
			res.push_back(temp);
		}
	}
    return res;
	
	}

void person::add_polymorphism(double r, double x, double y) {
	Gene temp1(r, x/2.0, y/2.0);
	Gene temp2(r, x/2.0, y/2.0);
	phen_x += x;
	phen_y += y;
	perf1 = calc_perf(phen_x, phen_y, 0.0, theta1);
	perf2 = calc_perf(phen_x, phen_y, archetype, theta2);
	gen1.push_back(temp1);
	gen2.push_back(temp2);

}

void person::set_niche(int n){
	niche = n;
}

void person::set_perf_1(double perf){
	perf1 = perf;
}

void person::set_perf_2(double perf){
	perf2 = perf;
}

void person::mutate(double mut_rate) {

	uniform_real_distribution<double> dist(0.0,1.0);
    uniform_real_distribution<double> angle(0.0,6.283185);
    //std::exponential_distribution<double> effectSize(1);
	std::gamma_distribution<double> effectSize(1, 1);
	uniform_int_distribution<int> DistZeroOne (0,1);
	poisson_distribution<int> mutNum(mut_rate);
	int mutation_num = mutNum(generator);


	for (int i = 0; i < mutation_num; i++) {
 
		double mut =  dist(generator);
		double ang = angle(generator);
		//double eff = 1;
		double eff = sqrt(effectSize(generator));
		//double choose = DistZeroOne(generator);
		//double eff = eff1*choose+eff2*(1-choose);

		Gene temp(mut, eff*cos(ang), eff*sin(ang));
		phen_x += temp.x;
		phen_y += temp.y;
		perf1 = calc_perf(phen_x, phen_y, 0.0, theta1);
		perf2 = calc_perf(phen_x, phen_y, archetype, theta2);
		if (DistZeroOne(generator) == 0) {
			gen1.push_back(temp);
		} else {
			gen2.push_back(temp);
		}
	}

}
pair <double,double> person::phenotype() {
    pair <double,double> ans(0.0,0.0);
    for (unsigned int i = 0 ; i < gen1.size(); i++) {
        ans.first += gen1.at(i).x;
        ans.second += gen1.at(i).y;
    }
    for (unsigned int i = 0 ; i < gen2.size(); i++) {
        ans.first +=gen2.at(i).x;
        ans.second +=gen2.at(i).y;
    }
	ans.first += common_pers.get_phen_x();
	ans.second += common_pers.get_phen_y();

    return ans;
}

string person::print() const{
    ostringstream a; 
    a << "gen1: \t";
    for (unsigned int i = 0; i <gen1.size(); i++){ 
        a << (gen1.at(i).r) << '\t';
        a << (gen1.at(i).x) << '\t';
        a << (gen1.at(i).y) << '\t';
    }
    a << "gen2: \t";
    for (unsigned int i = 0; i <gen2.size(); i++){
        a << (gen2.at(i).r) << '\t';
        a << (gen2.at(i).x) << '\t';
        a << (gen2.at(i).y) << '\t';
    }
    a << '\n';
    return a.str();
}

double person::get_phen_x() const{
	return phen_x;
}
double person::get_phen_y() const{
	return phen_y;
}

int person::get_niche(){
	return niche;
}
void person::set_phen_x(double x){
	phen_x = x;
}
void person::set_phen_y(double y){
	phen_y = y;
}
common_person::common_person(){
	phen_x = 0.0;
	phen_y = 0.0;
	pop_size = 0;
}

void common_person::reset_common_person(string filename){
	fixed.clear();
	phen_x = 0.0;
	phen_y = 0.0;
	ifstream infile_fixed;
	infile_fixed.open(filename);
	string line;
	std::getline(infile_fixed, line);

	std::vector<std::string> tokens;
	std::istringstream iss(line);
	std::string token;
	bool isgen1 = false;
	bool isgen2 = false;
	int order = 0;
	double r;
	double x;
	double y;
	while (std::getline(iss, token, '\t'))
		tokens.push_back(token);

	for (int i = 0; i < tokens.size(); i++){

	switch (order){
		case 0: //r
			r = stod(tokens.at(i));
			order += 1;
			break;
		case 1://x
			x = stod(tokens.at(i));
			order += 1;
			break;
		case 2://y
			y = stod(tokens.at(i));
			order = 0;
			Gene temp(r, x, y);
			fixed.insert(temp);
			phen_x += 2 * temp.x;
			phen_y += 2 * temp.y;
		}
	}
}
bool common_person::increase_gene(Gene g){
	if (common.count(g) == 0){
		common[g] = 0;
	}

	common[g] += 1;
	if (common[g] == 2 * pop_size){
		phen_x += 2*g.x;
		phen_y += 2*g.y;
		fixed.insert(g);
		return true;
	}
	return false;
}


int common_person::get_value(Gene g){
	return common[g];
}
void common_person::add_gene(Gene g){
	common[g] = 0;
}
void common_person::reset(){
	common.clear();
}
bool common_person::exist_element(Gene g){
	return (common.count(g)>0);
}
double common_person::get_phen_x(){
	return phen_x;
}
double common_person::get_phen_y(){
	return phen_y;
}

int common_person::get_pop_size(){
	return pop_size;
}

bool common_person::in_fixed(Gene g){
	return fixed.count(g) > 0;
}
void common_person::set_pop_size(const int popSize){
	pop_size = popSize;
}
string common_person::print() const{
	ostringstream a;
	unordered_set<Gene>::const_iterator it;
	for (it = fixed.begin(); it != fixed.end(); ++it)
	{
		a << (it->r) << '\t';
		a << (it->x) << '\t';
		a << (it->y) << '\t';
	}
	a << '\n';
	return a.str();
}