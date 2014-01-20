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

static profile 
constructProfile (const string &filename) {
    ifstream in(filename.c_str());
    if(!in)
	throw SMITHLABException("Could not open "+ filename);
    profile newProfile;
    newProfile.initialize();
    getline(in,newProfile.id);
    in>>newProfile.total_kmers;
    for(size_t j=0; j<smithlab::alphabet_size; ++j) {
	char atcg; // deal with ACGT
	in>>atcg;
	in>>newProfile.background[j];
    }
    size_t index, counts;
    while(in>>index) {
	in>>counts;
	newProfile.kmer_counts.insert(make_pair<size_t,size_t>(index,counts));
    }
    in.close();
    return newProfile;
}

//Calculate the probability of a word from it's decimal index
// static double 
// probWordFromIndex(size_t index, const vector<double> &background) {
//     double prob=1;
//     size_t base=smithlab::alphabet_size;
//     while(index != 0) {
// 	prob *= background[index%base];
// 	index /= base;
//     }
//     return prob;
// }
//
// static void
// computeInnerProduct(const vector<profile> &referenceVectors, 
// 	const profile &queryProfile,
// 	vector<double> &outFeature) {
// }
//
// static void
// computeD2shape(const vector<profile> &referenceVectors, 
// 	const profile &queryProfile,
// 	vector<double> &outFeature,
// 	const size_t k_value) {
// }

/**
 * Taking two profile, calculate their Jaccard similarity
 */
static double
computeJaccardIndex(const profile &referenceProfile,
	const profile &queryProfile) {

    size_t intersection = 0;
    size_t sizeOfUnion = 0;
    unordered_map<size_t,size_t> temp = queryProfile.kmer_counts;
    for(unordered_map<size_t, size_t>::const_iterator it(referenceProfile.kmer_counts.begin());
	    it != referenceProfile.kmer_counts.end(); ++it) {
	if(temp.find(it->first) != temp.end()) {
	    intersection += std::min(it->second,temp[it->first]);
	    temp[it->first] = std::max(it->second, temp[it->first]);
	}
	else
	    temp[it->first] = it->second;
    }
    for(unordered_map<size_t,size_t>::const_iterator it(temp.begin());
	    it != temp.end(); ++it) {
	sizeOfUnion += it->second;
    }
    return 1.0*intersection/sizeOfUnion;
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
	opt_parse.add_opt("distance",'d',"distance measure e.g. dotProduct, d2shape, jaccard",false,distance);

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
	
	vector<string> filenames;
	/* Read Bacteria Kmer Feature Vectors As Reference Files*/
	read_dir(dir,filename_suffix,filenames);

	for (size_t i = 0; i < input_filenames.size(); ++i) {
	    unordered_map<string,double> outFeature;
	    string queryID;
	    profile queryProfile = constructProfile(input_filenames[i]);
	    queryID = queryProfile.id;

	    for(size_t i=0; i < filenames.size(); ++i) {
		profile refProfile = constructProfile(filenames[i]);
		double val = 0.0;
		if(distance == "dotProduct"){ 
		    val = computeJaccardIndex(refProfile,queryProfile);
		}
		else if(distance == "d2shape"){ 
		    val = computeJaccardIndex(refProfile,queryProfile);
		}
		else if(distance == "jaccard"){ 
		    val = computeJaccardIndex(refProfile,queryProfile);
		}

		outFeature[refProfile.id] = val;
	    }

	    std::ofstream of;
	    if (!outfile.empty()) of.open(outfile.c_str(),
		    std::ofstream::out | std::ofstream::app);
	    std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
	    out << k_value << '\n';
	    out << queryID << '\n';
	    for (unordered_map<string,double>::const_iterator it(outFeature.begin());
		    it != outFeature.end(); ++it)
		out << it->first<<" "<< it->second<< '\n';
	}

	cout<<"The size of reference database is: " << filenames.size()<<endl;
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
