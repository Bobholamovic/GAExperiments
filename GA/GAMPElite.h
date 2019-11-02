#pragma once
#include "GAMP.h"

namespace GA
{
	namespace GAMPElite
	{
		// ��Ӣ�ุ���ӽ��ݻ��㷨
		class CEvolutionMPE: public GAMP::CEvolutionMP
		{
		protected:
			void SelectParents_(vector<CIndividual const*>&);
		public:
			int m_nCandidates = 2;	// ���ӿռ���ѡȡ�ĵ�����>=2
			int m_nElites = 2;	// ÿ�β����ӽ��ľ�Ӣ������
		public:
			CEvolutionMPE(
				const OBJ_FUNC_TYPE& objfunc, int n,
				int np, int nc,
				vector<double>&& lb, vector<double>&& ub,
				double pr, double pc
			) :
				CEvolutionMP(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc)
			{}

			CEvolutionMPE(
				const OBJ_FUNC_TYPE& objfunc, int n,
				int np, int nc,
				CHROM&& lb, CHROM&& ub,
				double pr, double pc,
				int m
			) :
				CEvolutionMP(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc, m)
			{}

			CEvolutionMPE(
				const OBJ_FUNC_TYPE& objfunc, int n,
				int np, int nc,
				CHROM&& lb, CHROM&& ub,
				double pr, double pc,
				int m, int k, int l
			) :
				CEvolutionMP(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc, m), 
				m_nElites(k), m_nCandidates(l)
			{}

			void Cross();
		};
	}
}