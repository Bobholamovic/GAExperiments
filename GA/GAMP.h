#pragma once
#include "GA.h"

namespace GA
{
	namespace GAMP
	{
		// 多父体杂交演化算法
		class CEvolutionMP : public NaiveGA::CEvolution
		{
		protected:
			CIndividual _Cross();
			virtual void _SelectParents(vector<CIndividual const*>&);
		public:
			int m_nParents = 4;	// 每次参与杂交的父本数
		public:
			CEvolutionMP(
				double(*objfunc)(CIndividual&), int n,
				int np, int nc,
				vector<double>&& lb, vector<double>&& ub,
				double pr, double pc
			) :
				CEvolution(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc)
			{}

			CEvolutionMP(
				double(*objfunc)(CIndividual&), int n,
				int np, int nc,
				vector<double>&& lb, vector<double>&& ub,
				double pr, double pc, 
				int m
			) :
				CEvolution(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc), m_nParents(m)
			{}

			void Cross();
		};
	}
}