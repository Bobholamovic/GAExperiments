#include "GAMPElite.h"

namespace GA
{
	namespace GAMPElite
	{
		void CEvolutionMPE::_SelectParents(vector<CIndividual const*>& parents)
		{
			// 选取 m_nElites 个最好的个体，另外 (m_nParents-m_nElites) 个个体从种群中随机选取
			for (auto i = 0; i < m_nParents; i++)
			{
				if (i < m_nElites)
				{
					// 由于 m_veciPopulation 有序，直接拷贝前m_nElites个
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
			// 收集候选人
			for (auto i = 0; i < m_nCandidates; i++)
			{
				candidates[i] = std::move(_Cross());
				// 决出前两名
				if (candidates[i].m_dFitValue >= dBestVal)	// 还是得加上等于的情况
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
			// 再次转移资源
			m_veciPopulation.push_back(std::move(candidates[nIdxBest]));
			m_veciPopulation.push_back(std::move(candidates[nIdxScnd]));
		}
	}
}