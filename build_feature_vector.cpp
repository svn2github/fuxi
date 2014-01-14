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

static void
computeInnerProduct(const vector<vector<double> > &referenceVectors, 
		const vector<double> &requestVector,
		vector<double> &outFeature) {

	for(vector<vector<double> >::const_iterator it = referenceVectors.begin();
			it != referenceVectors.end(); ++it) {
		double result=0;
		vector<double>::const_iterator it1(it->begin());
		vector<double>::const_iterator it2(requestVector.begin());
		while(it1 != it->end() && it2 != requestVector.end()) {
			result += (*it1)*(*it2);
			it1++;
			it2++;
		}
		if(it1 != it->end() || it2 != requestVector.end())
			throw SMITHLABException("The dimension of request vector "
					"and reference vector are not the same!");
		outFeature.push_back(result);
	}
}

static void
constructReferenceVector(const string &filename,
		const vector<vector<double> > &referenceVectors,
		vector<double> &outFeature,
		string &id) {

	vector<double> requestVector;
	ifstream in(filename.c_str());
	if(!in)
		throw SMITHLABException("Could not open "+ filename);
	string line;
	double kmer;
	getline(in,id);
	while(in>>kmer) {
		requestVector.push_back(kmer);
	}
	computeInnerProduct(referenceVectors,requestVector,outFeature);
}

int
main(int argc, const char **argv) {

	try {

		static const string filename_suffix = "cv";

		/****************** COMMAND LINE OPTIONS ********************/
		OptionParser opt_parse(strip_path(argv[0]),
				"makes feature vectors using inner product "
				"with reference genome from kmer vector file(s)",
				"<outfile> <dir> <infile1> [<infile2> ...]");

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
		vector<vector<double> > referenceVectors;
		vector<string> filenames;
		/* Read Bacteria Kmer Feature Vectors As Reference Files*/
		read_dir(dir,filename_suffix,filenames);

		/* Sort the reference genomes alphabetically, so that the reference can
		 * keep consistent under different systems
		 */
		sort(filenames.begin(),filenames.end());

		for(size_t i=0; i < filenames.size(); ++i) {
			vector<double> vectorForEachFile;
			ifstream in(filenames[i].c_str());
			if(!in)
				throw SMITHLABException("could not open file" + filenames[i]);
			string line;
			string id;
			double kmer;
			getline(in,id);
			while(in>>kmer)
				vectorForEachFile.push_back(kmer);
			referenceVectors.push_back(vectorForEachFile);
		}


		for (size_t i = 0; i < input_filenames.size(); ++i) {
			vector<double> outFeature;
			string requestID;
			constructReferenceVector(input_filenames[i],referenceVectors,
					outFeature,requestID);
			std::ofstream of;
			if (!outfile.empty()) of.open(outfile.c_str(),
					std::ofstream::out | std::ofstream::app);
			std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());

			out << requestID << '\n';
			for (vector<double>::iterator it(outFeature.begin());
					it != outFeature.end(); ++it)
				out << *it << '\n';
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
