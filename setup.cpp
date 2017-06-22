#include "setup.hpp"
#include "util.hpp"

void summary_stat::init(int dim)
{
     num_ld_scores=dim;
	 ld_score=d_createVector(dim);
}