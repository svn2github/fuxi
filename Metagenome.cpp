/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *                       Ehsan Behnam
 *
 *    Authors: Andrew D. Smith and Ehsan Behnam
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
#include "Metagenome.hpp"
#include "smithlab_utils.hpp"
#include "smithlab_os.hpp"

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <climits>
#include <tr1/unordered_map>

using std::string;
using std::vector;

typedef std::tr1::unordered_map<string, size_t> count_set;
typedef std::tr1::unordered_map<string, double> normalized_count_set;

//double
//Metagenome::get_value(const string &s) const {
//  count_set::const_iterator i(word_counts.find(s));
//  if (i == word_counts.end()) return 0.0;
//  else return i->second;
//}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////
///////   FUNCTIONS FOR I/O
//////
/////


string
Metagenome::tostring() const {
  std::ostringstream oss;
  oss << id << '\n';
  for (size_t i = 0; i < normalized_word_counts.size(); ++i)
    oss << normalized_word_counts[i] << '\n';
  return oss.str();
//  normalized_count_set::const_iterator i(normalized_word_counts.begin());
//  if (i != normalized_word_counts.end()) {
//    oss << i->first << '\t' << i->second;
//    ++i;
//    for (; i != normalized_word_counts.end(); ++i)
//      oss << '\n' << i->first << '\t' << i->second;
//  }
//  return oss.str();
}



/*IO functions are not finalized.
 * It is impossible now to reconstruct word_counts
 * from file since we just save normalized_word_counts.
 */
std::ostream&
operator<<(std::ostream &os, const Metagenome &mg) {
  return os << mg.tostring();
}


//std::istream&
//operator>>(std::istream &is, Metagenome &mg) {
//  is >> mg.id;
//
//  string tmp_word;
//  double tmp_count;
//  while (is >> tmp_word >> tmp_count)
//    mg.normalized_word_counts[tmp_word] = tmp_count;
//
//  return is;
//}
// inline static size_t



//static size_t
//get_kmer(const string::const_iterator a, string::const_iterator b) {
//   size_t multiplier = 1, index = 0;
//   do {
//     --b;
//     const size_t the_base = base2int(*b);
//     if (the_base == smithlab::alphabet_size)
//       return std::numeric_limits<size_t>::max();
//     index += base2int(*b)*multiplier;
//     multiplier *= smithlab::alphabet_size;
//   } while (b > a);
//   return index;
// }
//
//
//
//std::istream&
//operator>>(std::istream &is, Metagenome &mg) {
//  is >> mg.id;
//  mg.normalized_word_counts.clear();
//
//  //CAUTION: Assuming all 5-mers now!!!!!!!
//  for (size_t i = 0; i < pow(4.0, 5.0); ++i)
//    mg.normalized_word_counts.push_back(0.0);
//
//  string tmp_word;
//  double tmp_count;
//  while (is >> tmp_word >> tmp_count) {
//      size_t k = get_kmer(tmp_word.begin(), tmp_word.end());
//      if (k < 0 || k > 1023)
//        throw std::exception();
//      mg.normalized_word_counts[k] = tmp_count;
//  }
//
//  return is;
//}



std::istream&
operator>>(std::istream &is, Metagenome &mg) {
  is >> mg.id;
  mg.normalized_word_counts.clear();

  //CAUTION: Assuming all 5-mers now!!!!!!!
  for (size_t i = 0; i < pow(4.0, 5.0); ++i)
    mg.normalized_word_counts.push_back(0.0);

  double tmp_count;
  size_t k = 0;
  while (is >> tmp_count) {
      mg.normalized_word_counts[k] = tmp_count;
      ++k;
  }

  if (k != pow(4.0, 5.0))
    throw std::exception();

  return is;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////
///////   FUNCTIONS FOR COMPARING TWO METAGENOME OBJECTS
//////
/////



//double
//Metagenome::D2_distance(const Metagenome &other) const {
//  //basic D_2 for now. Note that it's a similarity not distance!
//  double score = 0.0;
//  for (count_set::const_iterator i(other.word_counts.begin());
//       i != other.word_counts.end(); ++i) {
//    count_set::const_iterator j(word_counts.find(i->first));
//    if (j != word_counts.end())
//      score += i->second*j->second;
//  }
//  return score;
//}



/* Function to compute the Euclidean distance between Metagenomes.
 */ 
//double
//Metagenome::euclidean_distance(const Metagenome &other) const {
//  double score = 0.0;
//  // first compute sum for coordinates non-zero in "other" vector
//  normalized_count_set::const_iterator i(other.normalized_word_counts.begin());
//  for (; i != other.normalized_word_counts.end(); ++i) {
//    normalized_count_set::const_iterator
//      j(normalized_word_counts.find(i->first));
//
//    const double diff = i->second -
//                      (j == normalized_word_counts.end() ? 0.0 : j->second);
//
//    score += diff*diff;
//  }
//  // then take care of coordinates that are zero on "other" vector
//  for (normalized_count_set::const_iterator i(normalized_word_counts.begin());
//       i != normalized_word_counts.end(); ++i) {
//
//      normalized_count_set::const_iterator
//        j(other.normalized_word_counts.find(i->first));
//
//    if (j == other.normalized_word_counts.end())
//      score += i->second*i->second;
//  }
//  return sqrt(score);
//}



//double
//norm(const normalized_count_set &mg) {
//  double ret = 0.0;
//  normalized_count_set::const_iterator i(mg.begin());
//  for (; i != mg.end(); ++i)
//    ret += i->second * i->second;
//  return sqrt(ret);
//}
static double
norm(const vector<double> &v) {
  double ret = 0.0;
  for (size_t i = 0; i < v.size(); ++i)
    ret += v[i] * v[i];
  return sqrt(ret);
}



double
Metagenome::compute_angle(const Metagenome &other) const {
  double angle = 0.0;

  assert(normalized_word_counts.size() == other.normalized_word_counts.size());
  for (size_t i = 0; i < normalized_word_counts.size(); ++i) {
      angle += normalized_word_counts[i] * other.normalized_word_counts[i];
  }

  angle /= (norm(normalized_word_counts) * norm(other.normalized_word_counts));

  // Under poor normalization rule, angle could be more than one if the points
  // are so close. This should be changed in the future.
  //
  if (abs(angle) > 1)
    angle = 1;
  //

  return acos(angle); // scale factor 180 / M_PI for "degree" conversion.
}
