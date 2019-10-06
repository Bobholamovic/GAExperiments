#include "GAMPElite.h"

namespace GA
{
	namespace GAMPElite
	{
		void CEvolutionMPE::_SelectParents(vector<CIndividual const*>& parents)
		{
			// ѡȡ m_nElites ����õĸ��壬���� (m_nParents-m_nElites) ���������Ⱥ�����ѡȡ
			for (auto i = 0; i < m_nParents; i++)
			{
				if (i < m_nElites)
				{
					// ���� m_veciPopulation ����ֱ�ӿ���ǰm_nElites��
					parents.push_back(&m_veciPopulation[i]);
				}
				else
				{
					parents.push_back(&m_veciPopulation[rand()%m_nPopSize]);
				}
			}

		}

		void CEvolutionMPE::Cross()
		{
			vector<CIndividual> candidates(m_nCandidates);
			int nIdxBest = 0, nIdxScnd = 0;
			double dBestVal = GA::MIN_VAL, dScndVal = GA::MIN_VAL;
			// �ռ���ѡ��
			for (auto i = 0; i < m_nCandidates; i++)
			{
				candidates[i] = std::move(_Cross());
				// ����ǰ����
				if (candidates[i].m_dFitValue >= dBestVal)	// ���ǵü��ϵ��ڵ����
				{
					nIdxScnd = nIdxBest;
					dScndVal = dBestVal;
					nIdxBest = i;
					dBestVal = candidates[i].m_dFitValue;
				}
				else if (candidates[i].m_dFitValue >= dScndVal)
				{
					nIdxScnd = i;
					dScndVal = candidates[i].m_dFitValue;
				}
			}
			// �ٴ�ת����Դ
			m_veciPopulation.push_back(std::move(candidates[nIdxBest]));
			m_veciPopulation.push_back(std::move(candidates[nIdxScnd]));
		}
	}
}