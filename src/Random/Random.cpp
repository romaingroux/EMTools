#include "Random.hpp"
#include <string>
#include <iostream>

bool rand_bernoulli(double p)
{   std::bernoulli_distribution dist(p) ;
    return dist(getRandomGenerator()) ;
}


std::vector<bool> rand_bernoulli(double p, size_t n)
{   std::vector<bool> vector(n) ;
    std::bernoulli_distribution dist(p) ;
    for(size_t i=0; i<n; i++)
    {   vector[i] = dist(getRandomGenerator()) ; }
    return vector ;
}


double rand_normal(double m, double sd)
{   std::normal_distribution<double> dist(m, sd) ;
    return dist(getRandomGenerator()) ;
}


std::vector<double> rand_normal(double m, double sd, double n)
{   std::vector<double> vector(n) ;
    std::normal_distribution<double> dist(m, sd) ;
    for(size_t i=0; i<n; i++)
    {   vector[i] = dist(getRandomGenerator()) ; }
    return vector ;
}

std::string rand_string(size_t length)
{    // the alphabet to sample
     static const std::string alphanum = "0123456789"
                                         "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                         "abcdefghijklmnopqrstuvwxyz" ;
     static const size_t lower = 0 ;
     static const size_t upper = alphanum.size() - 1 ;
     // sampling
     std::string s ;
     for(size_t i=0; i<length; i++)
     {   s.push_back(alphanum[rand_int_uniform(lower, upper)]) ; }
     return s ;
}
