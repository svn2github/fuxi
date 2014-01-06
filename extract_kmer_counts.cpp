/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *                       Ehsan Behnam
 *                       Wenzheng Li
 *
 *    Authors: Andrew D. Smith, Ehsan Behnam and Wenzheng Li
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
#include <iostream>
#include <tr1/unordered_map>
#include <cmath>

#include "OptionParser.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"
#include "Metagenome.hpp"

using std::string;
using std::vector;
using std::cerr;
using std::cout;
using std::endl;
using std::tr1::unordered_map;



struct stats_background {
   vector <double> background;
   size_t n_characters;
   stats_background() {}
   void initialize();
};

void
stats_background::initialize() {
   n_characters = 0;
   for (size_t i = 0; i < smithlab::alphabet_size; ++i)
      background.push_back(0.0);
}

struct stats_collector {
   stats_collector() : bad_kmers(0), total_kmers(0) {}
   size_t bad_kmers;
   size_t total_kmers;
};

std::ostream&
operator<<(std::ostream &os, const stats_collector &sc) {
   return os << sc.total_kmers << ", " << sc.bad_kmers;
}

inline static size_t
get_kmer(const string::const_iterator a, string::const_iterator b) {
   size_t multiplier = 1, index = 0;
   do {
      --b;
      const size_t the_base = base2int(*b);
      if (the_base == smithlab::alphabet_size)
         return std::numeric_limits<size_t>::max();
      index += base2int(*b)*multiplier;
      multiplier *= smithlab::alphabet_size;
   } while (b > a);
   return index;
}

static void
extract_kmer_counts_sequence(const int k_value, const string &sequence,
      unordered_map<size_t, size_t> &kmer_counts,
      stats_collector &sc) {
   const string::const_iterator limit = sequence.end() - k_value + 1;
   for (string::const_iterator i(sequence.begin()); i != limit; ++i) {
      const bool good_kmer = (count(i, i + k_value, 'A')+count(i, i + k_value, 'G')+count(i, i + k_value, 'C')+
                  count(i, i + k_value, 'T')+count(i, i + k_value, 'a')+count(i, i + k_value, 'g')
                  +count(i, i + k_value, 'c')+count(i, i + k_value, 't')==k_value);
      if (good_kmer) {
         const size_t the_kmer=get_kmer(i,i+k_value); //get_kmer(i, i + k_value);
         unordered_map<size_t, size_t>::iterator j = kmer_counts.find(the_kmer);
         if (j == kmer_counts.end())
            kmer_counts[the_kmer] = 1;
         else
            j->second++;
      }
      else ++sc.bad_kmers;
      ++sc.total_kmers;
   }
}

static void
update_background(string &line, stats_background &b) {
   for (size_t i = 0; i < line.length(); ++i)
      switch (line[i]) {
         case 'A':
            ++b.background[0];
            ++b.n_characters;
            break;

         case 'C':
            ++b.background[1];
            ++b.n_characters;
            break;

         case 'G':
            ++b.background[2];
            ++b.n_characters;
            break;

         case 'T':
            ++b.background[3];
            ++b.n_characters;
            break;
            //

         case 'a':
            ++b.background[0];
            ++b.n_characters;
            break;

         case 'c':
            ++b.background[1];
            ++b.n_characters;
            break;
         case 'g':
            ++b.background[2];
            ++b.n_characters;
            break;

         case 't':
            ++b.background[3];
            ++b.n_characters;
            break;

      }
}

static double
calc_prob_word(string &word, const stats_background &b) {
   double ret = 1.0;
   for (size_t i = 0; i < word.length(); ++i)
      switch (word[i]) {
         case 'A': ret *= b.background[0]; break;
         case 'C': ret *= b.background[1]; break;
         case 'G': ret *= b.background[2]; break;
         case 'T': ret *= b.background[3]; break;
                   //
         case 'a': ret *= b.background[0]; break;
         case 'c': ret *= b.background[1]; break;
         case 'g': ret *= b.background[2]; break;
         case 't': ret *= b.background[3]; break;
      }
   return ret;
}


static bool
is_fastq_sequence_line(const size_t line_number) {
  return (line_number % 4 == 1);
}

static void
extract_kmer_counts_from_fastq(const int k_value,
                               const vector<string> &input_filenames,
                               unordered_map<size_t, size_t> &kmer_counts,
                               stats_collector &sc, size_t &read_length,
                               vector<double> &frequency) {

   stats_background b;
   b.initialize();
   read_length = 0;

   for (size_t i = 0; i < input_filenames.size(); i++) {
      std::ifstream in(input_filenames[i].c_str());
      if (!in)
         throw SMITHLABException("could not open " + input_filenames[i]);

      string line;
      size_t line_count = 0;
      while (getline(in, line)) {
         if (is_fastq_sequence_line(line_count)) {
             extract_kmer_counts_sequence(k_value, line, kmer_counts, sc);
             update_background(line, b);
             read_length += (line.length()-k_value+1);
         }
         ++line_count;
      }
   }

   //used for debug
   cout<<"read_length:"<<read_length<<endl;
   //static const size_t d = pow(smithlab::alphabet_size, k_value);
   //const double n_reads = read_length;

   for (size_t i = 0; i < smithlab::alphabet_size; ++i) {
      b.background[i] /= b.n_characters;
      frequency.push_back(b.background[i]);
   }
   //unordered_map<string, size_t>::const_iterator i = kmer_counts.begin();

   cout<<"bad_kmers:"<<sc.bad_kmers<<endl;
}



static void
extract_kmer_counts_from_fasta(const int k_value,
                               const vector<string> &input_filenames,
                               unordered_map<size_t, size_t> &kmer_counts,
                               stats_collector &sc,
                               size_t &read_length,
                               vector<double> &frequency) {

   stats_background b;
   b.initialize();
   read_length = 0;

   for (size_t i = 0; i < input_filenames.size(); i++) {
      std::ifstream in(input_filenames[i].c_str());
      if (!in)
         throw SMITHLABException("could not open " + input_filenames[i]);

       string description;
	   getline(in,description);
	   string line;
	   string seq;
	   while (getline(in, line)) {
	     seq+=line;
	   }
	   extract_kmer_counts_sequence(k_value,seq, kmer_counts, sc);
	   update_background(seq, b);
	   read_length += (seq.length()-k_value+1);
   }

   //used for debug
   cout<<"read_length:"<<read_length<<endl;
   static const size_t d = pow(smithlab::alphabet_size, k_value);
   const double n_reads = read_length;

   for (size_t i = 0; i < smithlab::alphabet_size; ++i) {
      b.background[i] /= b.n_characters;
      frequency.push_back(b.background[i]);
   }
   //unordered_map<string, size_t>::const_iterator i = kmer_counts.begin();

   cout<<"bad_kmers:"<<sc.bad_kmers<<endl;
}



int
main(int argc, const char **argv) {

   try {

      bool VERBOSE = false, isFASTA = false;
      int k_value = 0;
      string id;

      /*********************** COMMAND LINE OPTIONS **************************/
      OptionParser opt_parse(strip_path(argv[0]),
            "makes k-mer counts file from FASTQ/FASTA file(s)",
            "<outfile> <infile1> [<infile2> ...]");
      opt_parse.add_opt("id", 'i', "metagenome id", true, id);
      opt_parse.add_opt("kmer", 'k', "word size", true, k_value);
      opt_parse.add_opt("fasta_flag", 'f', "fasta or fastq", false, isFASTA);
      opt_parse.add_opt("verbose", 'v', "print more run info", false, VERBOSE);

      vector<string> leftover_args;
      opt_parse.parse(argc, argv, leftover_args);
      if (argc < 2 || opt_parse.help_requested()) {
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
      if (leftover_args.size() < 2) {
         cerr << opt_parse.help_message() << endl;
         return EXIT_SUCCESS;
      }
      const string outfile(leftover_args.front());
      vector<string> input_filenames;
      copy(leftover_args.begin() + 1, leftover_args.end(),
            back_inserter(input_filenames));
      /*********************** END COMMAND LINE OPTIONS **********************/

      if (VERBOSE) {
         std::clog << "input files:" << endl;
         copy(input_filenames.begin(), input_filenames.end(),
               std::ostream_iterator<string>(std::clog, "\n"));
      }

      unordered_map<size_t, size_t> kmer_counts;
      stats_collector sc;

      size_t total_length=0;
      std::vector<double> background;
      if (isFASTA)
        extract_kmer_counts_from_fasta(k_value, input_filenames,
                                        kmer_counts, sc, total_length,
                                        background);
      else
        extract_kmer_counts_from_fastq(k_value, input_filenames,
                                        kmer_counts, sc, total_length,
                                        background);

      /************** WRITING THE k-CV WITHOUT CENTRALIZATION ****************/
      std::ofstream of;
      if (!outfile.empty()) of.open(outfile.c_str());
      std::ostream out(outfile.empty() ? std::cout.rdbuf() : of.rdbuf());
      out << id << endl;
      out << total_length-sc.bad_kmers << endl;
      out << "A " << background[0] << endl;
      out << "C " << background[1] << endl;
      out << "G " << background[2] << endl;
      out << "T " << background[3] << endl;

      map<size_t,size_t> ordered_kmers(kmer_counts.begin(),kmer_counts.end());
      map<size_t,size_t>::const_iterator it = ordered_kmers.begin();
      while(it != ordered_kmers.end()) {
         out<<it->first<<" "<<it->second<<'\n';
         ++it;
      }

      if (VERBOSE)
         std::clog << "(total kmers, bad kmers) = " << sc << endl;
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
