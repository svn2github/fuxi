/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2013 University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Ehsan Behnam, Andrew D. Smith
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
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef METAGENOME_HPP
#define METAGENOME_HPP

#include <string>
#include <vector>
#include <numeric>
#include <tr1/unordered_map>

class Metagenome {
public:
  /* Copy constructor
   * assignment operator
   * tostring() const
   */ 
  Metagenome() {}
  Metagenome(const std::string &id_in,
	     const std::tr1::unordered_map<std::string, size_t>
             &word_counts_in,
             const std::vector<double>
             &normalized_word_counts_in) :
    id(id_in), word_counts(word_counts_in),
    normalized_word_counts(normalized_word_counts_in) {}
  
  // Accessors
  bool good() const {return !id.empty() && !normalized_word_counts.empty();}
  std::string get_id() const {return id;}
  size_t size() const {return normalized_word_counts.size();}
  std::string tostring() const;
  std::vector<double>::const_iterator begin() const
    {return normalized_word_counts.begin();}

  inline double operator[](size_t i) const {return normalized_word_counts[i];}

  //double get_value(const std::string &) const;
  
  // Functions for comparing two Metagenomes
  //double D2_distance(const Metagenome &other) const;
  //double euclidean_distance(const Metagenome &other) const;
  double compute_angle(const Metagenome &other) const;
  double compute_dot_product(const std::vector<double> &v) const
    {return std::inner_product(normalized_word_counts.begin(),
                               normalized_word_counts.end(),
                               v.begin(), 0.0);}
  friend std::istream&
  operator>>(std::istream &, Metagenome &);
  
private:
  std::string id;
  std::tr1::unordered_map<std::string, size_t> word_counts;
  //std::tr1::unordered_map<std::string, double> normalized_word_counts;
  std::vector<double> normalized_word_counts;
};

std::ostream&
operator<<(std::ostream &, const Metagenome &);

#endif
