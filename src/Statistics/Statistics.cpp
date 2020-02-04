#include <Statistics.hpp>
#include <cmath>   // M_PI, pow(), sqrt(), log(), lgamma()
#include <boost/math/distributions.hpp>  // beta_distribution
#include <boost/random.hpp>


double normal_pmf(double x, double mean, double sd)
{   static double pi_2 = 2.*M_PI ;
    return ( 1. / ( sd * sqrt(pi_2) )) * exp(-0.5 * pow((x-mean)/sd, 2.0 ) );
}

double poisson_pmf(int x, double lambda)
{   if((x != 0) and (lambda == 0))
    {   return 0. ; }
    else if((x == 0) and (lambda == 0))
    {   return 1. ; }

    // if(x + lambda == 0)
    // {   return 1. ; }
    return exp(x * log(lambda) - lgamma(x + 1.0) - lambda);
}

double poisson_cdf(int x, double lambda, bool lower_tail)
{   double prob = 1. ;
    if(lambda != 0)
    {   boost::math::poisson_distribution<> pois_dist(lambda) ;
        prob = cdf(pois_dist, x) ;
    }
    // prob of an observation as big as x
    if(lower_tail)
    {   return prob ; }
    // prob of an observation at least as big as x (pvalue)
    else
    {   return 1 - prob ; }
}

double beta_pmf(double x, double alpha, double beta)
{
    boost::math::beta_distribution<> beta_dist(alpha, beta) ;
    double y = quantile(beta_dist, x) ;

    return y ;
}


