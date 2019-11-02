#pragma once
#include "GAMPElite.h"

namespace GA
{
	namespace GLME
	{
		// ȫ��-�ֲ��ݻ��㷨
		class CGLME : public GAMP::CEvolutionMP
		{
		private:
			int m_nGlobalGen__ = 10;	// ȫ���ݻ�����
			int m_nCentersMax__ = 10;	// ���ظ���
			int m_nCentersMin__ = 10;	// ��С�ظ���
			double m_dEpsilon = 0.5;	// ����ѡ�������epsilon
			double m_dAlpha = 1.0;		// ����ѡ��������ӣ�alpha
			int m_nSubPopSize__ = 50;	// ÿ��������Ⱥ��С
			vector<GAMPElite::CEvolutionMPE> m_veciClusters__;	// �صļ���

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
			int m_nCandidates = 2;	// ���ӿռ���ѡȡ�ĵ�����>=2
			int m_nElites = 2;	// ÿ�β����ӽ��ľ�Ӣ������

		private:
			static double CalcDist__(const CHROM&, const CHROM&);	// ��������Ⱦɫ��֮��ľ���
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

			// ʵ���ԣ������ݻ�
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
