/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *                       Ehsan Behnam
 *                       Wenzheng Li
 *
 *    Authors: Andrew D. Smith , Ehsan Behnam and Wenzheng Li
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 */

#include <string>
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <tr1/unordered_map>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Metagenome.hpp"


using std::string;
using std::ifstream;
using std::vector;
using std::cerr;
using std::cout;
using std::cin;
using std::endl;
using std::flush;
using std::make_pair;
using std::tr1::unordered_map;

struct profile {
    string id; //The ID of the genome
    size_t total_kmers; // The number of total kmers;
    vector <double> background; // The frequency of A,C,G,T respectively
    unordered_map<size_t, size_t>  kmer_counts; // The pairs of kmers
    void initialize();
    void print();
};

void
profile::initialize() {
    total_kmers = 0;
    for (size_t i=0; i<smithlab::alphabet_size; ++i) 
	background.push_back(0.0);
}

void
profile::print() {
    cout<<"The profile for **"<< id <<"**:\n";
    cout<<"total_kmers:"<<total_kmers<<"\n";
    cout<<"A:"<<background[0]<<"\n";
    cout<<"C:"<<background[1]<<"\n";
    cout<<"G:"<<background[2]<<"\n";
    cout<<"T:"<<background[3]<<"\n";
    cout<<"Kmer_counts:\n";
    unordered_map<size_t,size_t>::iterator it = kmer_counts.begin();
    while(it != kmer_counts.end()) {
	cout<<it->first<<" "<<it->second<<'\n';
	++it;
    }
}

//Calculate the probability of a word from it's decimal index
static double 
probWordFromIndex(size_t index, const vector<double> &background) {
    double prob=1;
    size_t base=smithlab::alphabet_size;
    while(index != 0) {
	prob *= background[index%base];
	index /= base;
    }
    return prob;
}

static void
computeInnerProduct(const vector<profile> &referenceVectors, 
	const profile &requestVector,
	vector<double> &outFeature) {
}

static void
computeD2shape(const vector<profile> &referenceVectors, 
	const profile &requestProfile,
	vector<double> &outFeature,
	const size_t k_value) {

    //static const size_t dim = pow(smithlab::alphabet_size,k_value);
    //Processing request kmer vector by reducing mean value
    //for(size_t i=0; i<dim; ++i) {
    //    double prob=probWordFromIndex(i, requestProfile.background);
    //    unordered_map<size_t,size_t>::iterator it = (requestProfile.kmer_counts).find(i);
    //    if(it == requestProfile.kmer_counts.end())
    //        requestProfile.kmer_counts[i] = 0;
    //    requestProfile.kmer_counts[i] -= (requestProfile.total_kmers*prob);
    //}

    //for(vector<profile>::const_iterator it = referenceVectors.begin();
    //        it != referenceVectors.end(); ++it) {
    //    //Processing reference bacteria kmer vector by reducing mean value
    //    for(size_t i=0; i<dim; ++i) {
    //        double prob=probWordFromIndex(i, it->background);
    //        unordered_map<size_t,size_t>::iterator it_n = it->kmer_counts.find(i);
    //        if(it_n == it->kmer_counts.end())
    //    	it->kmer_counts[i] = 0;
    //        it->kmer_counts[i] -= (it->total_kmers*prob);
    //    }
    //}
}

static void
computeJaccardIndex(const vector<profile> &referenceVectors, 
	const profile &requestVector,
	vector<double> &outFeature) {
}

static void
constructReferenceVector(const string &filename,
	const vector<profile> &referenceVectors,
	vector<double> &outFeature,
	string &id, const size_t k_value,const string &distance) {
    ifstream in(filename.c_str());
    if(!in)
	throw SMITHLABException("Could not open "+ filename);
    profile requestProfile;
    requestProfile.initialize();
    getline(in,requestProfile.id);
    in>>requestProfile.total_kmers;
    for(size_t j=0; j<smithlab::alphabet_size; ++j) {
	char atcg; // deal with ACGT
	in>>atcg;
	in>>requestProfile.background[j];
    }
    size_t index, counts;
    while(in>>index) {
	in>>counts;
	requestProfile.kmer_counts.insert(make_pair<size_t,size_t>(index,counts));
    }
    if(distance == "dotProduct")
	computeInnerProduct(referenceVectors,requestProfile,outFeature);
    else if(distance == "d2shape")
	computeD2shape(referenceVectors,requestProfile,outFeature,k_value);
    else if(distance == "jaccard") 
	computeJaccardIndex(referenceVectors,requestProfile,outFeature);
}

int
main(int argc, const char **argv) {

    try {

	static const string filename_suffix = "cv";
	size_t k_value = 0;
	string distance = "dotProduct";

	/****************** COMMAND LINE OPTIONS ********************/
	OptionParser opt_parse(strip_path(argv[0]),
		"makes feature vectors using inner product "
		"with reference genome from kmer vector file(s)",
		"<outfile> <dir> <infile1> [<infile2> ...]");
	opt_parse.add_opt("kmer",'k',"word size",true,k_value);
	opt_parse.add_opt("distance",'d',"distance measure",true,distance);

	vector<string> leftover_args;
	opt_parse.parse(argc, argv, leftover_args);
	if (argc == 1 || opt_parse.help_requested()) {
	    cerr << opt_parse.help_message() << endl
		<< opt_parse.about_message() << endl;
	    return EXIT_SUCCESS;
	}
	if (opt_parse.about_requested()) {
	    cerr << opt_parse.about_message() << endl;
	    return EXIT_SUCCESS;
	}
	if (opt_parse.option_missing()) {
	    cerr << opt_parse.option_missing_message() << endl;
	    return EXIT_SUCCESS;
	}
	if (leftover_args.size() < 3) {
	    cerr << opt_parse.help_message() << endl;
	    return EXIT_SUCCESS;
	}
	const string outfile(leftover_args.front());
	const string dir(leftover_args[1]);
	vector<string> input_filenames;
	copy(leftover_args.begin() + 2, leftover_args.end(),
		back_inserter(input_filenames));
	/****************** END COMMAND LINE OPTIONS *****************/

	/////////////////////////////////////////////////////////////////
	// READING IN THE DATABASE OF BACTERIAL GENOME FEATURES AND THEY ARE
	// USED TO CONSTRUCT NEW FEATURES FOR OTHER METAGENOME
	vector<profile> referenceVectors;
	vector<string> filenames;
	/* Read Bacteria Kmer Feature Vectors As Reference Files*/
	read_dir(dir,filename_suffix,filenames);

	/* Sort the reference genomes alphabetically, so that the reference can
	 * keep consistent under different systems
	 */
	sort(filenames.begin(),filenames.end());


	// If k_value is small, we can load all of the reference bacteria kmer vectors, the
	// memory needed here is about 1 gigabytes for k=6
	if(k_value <= 6) {
	    for(size_t i=0; i < filenames.size(); ++i) {
		profile eachProfile;
		eachProfile.initialize();

		ifstream in(filenames[i].c_str());
		if(!in)
		    throw SMITHLABException("could not open file" + filenames[i]);
		getline(in,eachProfile.id);
		in>>eachProfile.total_kmers;
		for(size_t j=0; j<smithlab::alphabet_size; ++j) {
		    char atcg; // deal with ACGT
		    in>>atcg;
		    in>>eachProfile.background[j];
		}
		size_t index, counts;
		while(in>>index) {
		    in>>counts;
		    eachProfile.kmer_counts.insert(make_pair<size_t,size_t>(index,counts));
		}

		//eachProfile.print();
		referenceVectors.push_back(eachProfile);
	    }


	    for (size_t i = 0; i < input_filenames.size(); ++i) {
		vector<double> outFeature;
		string requestID;
		constructReferenceVector(input_filenames[i],referenceVectors,
			outFeature,requestID,k_value,distance);
		std::ofstream of;
		if (!outfile.empty()) of.open(outfile.c_str(),
			std::ofstream::out | std::ofstream::app);
		std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

		out << k_value << '\n';
		out << requestID << '\n';
		for (vector<double>::iterator it(outFeature.begin());
			it != outFeature.end(); ++it)
		    out << *it << '\n';
	    }
	}
    }
    catch (const SMITHLABException &e) {
	cerr << e.what() << endl;
	return EXIT_FAILURE;
    }
    catch (std::bad_alloc &ba) {
	cerr << "ERROR: could not allocate memory" << endl;
	return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
