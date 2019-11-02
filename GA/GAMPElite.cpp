#include "GAMPElite.h"

namespace GA
{
	namespace GAMPElite
	{
		void CEvolutionMPE::SelectParents_(vector<CIndividual const*>& parents)
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
			//// FIXME: 此法似乎很不安全，容易move同一资源两次，考虑换成排序吧
			//int nIdxBest = 0, nIdxScnd = 1;
			//double dBestVal = GA::MIN_VAL-1, dScndVal = GA::MIN_VAL-2;	// 加上一个偏移值以应对MIN_VAL情况
			//// 收集候选人
			//for (auto i = 0; i < m_nCandidates; i++)
			//{
			//	candidates[i] = std::move(_Cross());
			//	// 决出前两名
			//	if (candidates[i].m_dFitValue >= dBestVal)	// 还是得加上等于的情况
			//	{
			//		nIdxScnd = nIdxBest;
			//		dScndVal = dBestVal;
			//		nIdxBest = i;
			//		dBestVal = candidates[i].m_dFitValue;
			//	}
			//	else if (candidates[i].m_dFitValue >= dScndVal)
			//	{
			//		nIdxScnd = i;
			//		dScndVal = candidates[i].m_dFitValue;
			//	}
			//}
			//// 再次转移资源
			//m_veciPopulation.push_back(std::move(candidates[nIdxBest]));
			//m_veciPopulation.push_back(std::move(candidates[nIdxScnd]));
			for (auto &e: candidates)
				e = std::move(Cross_());
			sort(candidates.begin(), candidates.end(), this->Criterion_);
			m_veciPopulation.push_back(std::move(candidates[0]));
			m_veciPopulation.push_back(std::move(candidates[1]));
		}
	}
}