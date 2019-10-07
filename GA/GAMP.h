#pragma once
#include "GA.h"

namespace GA
{
	namespace GAMP
	{
		// �ุ���ӽ��ݻ��㷨
		class CEvolutionMP : public NaiveGA::CEvolution
		{
		protected:
			CIndividual _Cross();
			virtual void _SelectParents(vector<CIndividual const*>&);
		public:
			int m_nParents = 4;	// ÿ�β����ӽ��ĸ�����
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