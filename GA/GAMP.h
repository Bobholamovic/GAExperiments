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
			CIndividual Cross_();
			virtual void SelectParents_(vector<CIndividual const*>&);
		public:
			int m_nParents = 4;	// ÿ�β����ӽ��ĸ�����
		public:
			CEvolutionMP(
				const OBJ_FUNC_TYPE& objfunc, int n,
				int np, int nc,
				CHROM&& lb, CHROM&& ub,
				double pr, double pc
			) :
				CEvolution(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc)
			{}

			CEvolutionMP(
				const OBJ_FUNC_TYPE& objfunc, int n,
				int np, int nc,
				CHROM&& lb, CHROM&& ub,
				double pr, double pc, 
				int m
			) :
				CEvolution(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc), m_nParents(m)
			{}

			void Cross();
		};
	}
}