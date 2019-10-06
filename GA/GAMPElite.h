#pragma once
#include "GAMP.h"

namespace GA
{
	namespace GAMPElite
	{
		// 多父体杂交演化算法
		class CEvolutionMPE: public GAMP::CEvolutionMP
		{
		protected:
			void _SelectParents(vector<CIndividual const*>&);
		public:
			int m_nCandidates = 2;	// 在子空间中选取的点数，>=2
			int m_nElites = 2;	// 每次参与杂交的精英父本数
		public:
			CEvolutionMPE(
				double(*objfunc)(CIndividual&), int n,
				int np, int nc,
				vector<double>&& lb, vector<double>&& ub,
				double pr, double pc
			) :
				CEvolutionMP(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc)
			{}

			CEvolutionMPE(
				double(*objfunc)(CIndividual&), int n,
				int np, int nc,
				vector<double>&& lb, vector<double>&& ub,
				double pr, double pc,
				double m
			) :
				CEvolutionMP(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc, m)
			{}

			CEvolutionMPE(
				double(*objfunc)(CIndividual&), int n,
				int np, int nc,
				vector<double>&& lb, vector<double>&& ub,
				double pr, double pc,
				double m, double k, double l
			) :
				CEvolutionMP(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc, m), 
				m_nElites(k), m_nCandidates(l)
			{}

			void Cross();
		};
	}
}