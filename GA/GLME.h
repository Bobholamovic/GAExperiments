#pragma once
#include "GAMPElite.h"

namespace GA
{
	namespace GLME
	{
		// 全局-局部演化算法
		class CGLME : public GAMP::CEvolutionMP
		{
		private:
			int m_nGlobalGen__ = 10;	// 全局演化代数
			int m_nCentersMax__ = 10;	// 最大簇个数
			int m_nCentersMin__ = 10;	// 最小簇个数
			double m_dEpsilon = 0.5;	// 中心选择参数，epsilon
			double m_dAlpha = 1.0;		// 中心选择比例因子，alpha
			int m_nSubPopSize__ = 50;	// 每个簇中种群大小
			vector<GAMPElite::CEvolutionMPE> m_veciClusters__;	// 簇的集合

		private:
			class Constraint__
			{
			public:
				double radius = 0.0;
				CHROM center;
			public:
				Constraint__(double r, const CHROM& c):
					radius(r), center(c)
				{}
				double operator()(CHROM& ind)
				{
					// ||X-Xi|| <= radius
					return radius - CGLME::CalcDist__(ind, this->center);
				}
			};

		public:
			int m_nCandidates = 2;	// 在子空间中选取的点数，>=2
			int m_nElites = 2;	// 每次参与杂交的精英父本数

		private:
			static double CalcDist__(const CHROM&, const CHROM&);	// 计算两个染色体之间的距离
			static vector<CIndividual> GetAllCenters__(const vector<CIndividual>&, double);
			using MPE_ITR = vector<GAMPElite::CEvolutionMPE>::iterator;
			static void LocEvolve(MPE_ITR, MPE_ITR, int);

		public:
			CGLME(
				const OBJ_FUNC_TYPE& objfunc, int n,
				int np, int nc,
				vector<double>&& lb, vector<double>&& ub,
				double pr, double pc,
				int m, int k, int l, 
				int gg
			) :
				CEvolutionMP(objfunc, n, np, nc, std::move(lb), std::move(ub), pr, pc, m),
				m_nElites(k), m_nCandidates(l), m_nGlobalGen__(gg)
			{}

			void FormClusters();
			vector<CIndividual> PickSolutions(double);
			void NextGeneration();

			// 实验性，并行演化
			void ParEvolve(int);

			inline int GetGlobalGen() const
			{
				return m_nGlobalGen__;
			}

			inline int GetSubPopSize() const
			{
				return m_nSubPopSize__;
			}

			inline void ConfigLocalParams(int min, int max, double eps, double alpha, int n1)
			{
				m_nCentersMax__ = max;
				m_nCentersMin__ = min;
				m_dEpsilon = eps;
				m_dAlpha = alpha;
				m_nSubPopSize__ = n1;
			}
		};
	}
}
