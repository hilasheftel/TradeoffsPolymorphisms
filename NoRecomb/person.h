#pragma once

#include <stdlib.h>
#include <vector>
#include <random>
#include <sstream>
#include <string>
#include <math.h>  
#include <unordered_map>
#include <unordered_set>


using namespace std;

extern default_random_engine generator;
extern double archetype;
extern double theta;
extern double lambda;
extern double PI;

class Gene{
public:
	double r;
	double x;
	double y;

	Gene(double r, double x, double y) : r(r), x(x), y(y) {};
	bool operator==(const Gene& other) const{
		return ((r == other.r) &&
			(x == other.x) &&
			(y == other.y));
	};
};

namespace std{
	template<>
	struct hash<Gene>{
		size_t operator()(Gene const& g) const{
			return hash<double>()(g.r);
		}

	};
}

typedef vector<Gene> Chromosome;

/*
This class contains the information about the individuals reproductive abilities and mutations.
*/

class common_person
{
private:
	unordered_map<Gene, int> common;
	unordered_set<Gene> fixed;
	double phen_x;
	double phen_y;
	int pop_size;
public:
	common_person();
	void reset_common_person(string);
	bool increase_gene(Gene);
	int get_value(Gene);
	void set_pop_size(const int);
	int get_pop_size();
	void add_gene(Gene);
	void add_to_fixed(Gene);
	void reset();
	bool exist_element(Gene);
	double get_phen_x();
	double get_phen_y();
	string print() const;
	bool in_fixed(Gene);

};

class person {
	//gen1 and gen2 are the 2 copies of the genome (diploids) and their effect on the phenotype.
private:
	int niche;
	Chromosome gen1;
	Chromosome gen2;
	double phen_x;
	double phen_y;
	//For reccombination, use after uniting 2 genes for infinite recombination
	struct {
		bool operator()(Gene& t1, Gene& t2)
		{
			return (t1.r < t2.r);
		}
	} sortGenome;

	//bool sortGenome (Gene t1,Gene t2) { return (t1.at(0)<t2.at(0)); }

public:
	double perf1;
	double perf2;
	//born as wiltype
	person(int);
	//born as recombination
	person(int, Chromosome, Chromosome);
	person(string);
	// mutation inducement
	//void mutate(double, int);
	//preparing the gametes
	//preparing the gametes using infinitie recommbination
	Chromosome ForRecombination(double);
	Chromosome choosecopy(double);

	Chromosome ForInfiniteRecombination(double);
	//mapping the genotype to phenotype.
	pair <double, double> phenotype();
	void mutate(double);
	void add_polymorphism(double, double, double);
	//print person
	string print() const;
	Chromosome get_gene(int);
	void set_niche(int);
	void set_perf_1(double);
	void set_perf_2(double);
	void set_phen_x(double);
	void set_phen_y(double);
	double get_phen_x() const;
	double get_phen_y() const;
	int get_niche();

	bool operator==(const person& other) const{
		return ((phen_x == other.phen_x) &&
			(phen_y == other.phen_y));
	};
};

size_t hash_combine(size_t lhs, size_t rhs);

namespace std{
	template<>
	struct hash<person>{
		size_t operator()(person const& p) const{
			return hash_combine(hash<double>()(p.get_phen_x()), hash<double>()(p.get_phen_y()));
			// return hash<pair<double,double>>()(pair<double, double> (p.perf1, p.perf2));
		}

	};

	template<>
	struct hash<pair<int, int>>{
		size_t operator()(pair<int, int> p) const{
			return hash_combine(hash<int>()(p.first), hash<int>()(p.second));
			// return hash<pair<double,double>>()(pair<double, double> (p.perf1, p.perf2));
		}

	};

}


