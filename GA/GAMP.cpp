#include "GAMP.h"

namespace GA
{
	namespace GAMP
	{
		void CEvolutionMP::SelectParents_(vector<CIndividual const*>& parents)
		{
			// 等概率随机选取 m_nParents 个父本
			for (auto i = 0; i < m_nParents; i++)
			{
				parents.push_back(&m_veciPopulation[rand()%m_nPopSize]);
			}
		}

		CIndividual CEvolutionMP::Cross_()
		{
			// 多父体杂交
			vector<double> vecdAlpha(m_nParents,0.0);
			CIndividual child;
			double dLb = -0.5, dUb = 1.5;	// 选取alpha的上下界
			double dSum = 0.0;
			vector<CIndividual const *> veciParents;

			SelectParents_(veciParents);	// 挑选父本

			// 随机设置 alpha
			for (auto i = 0; i < m_nParents-1; i++)
			{
				vecdAlpha[i] = RandRange((dLb < -0.5 ? -0.5 : dLb), (dUb > 1.5 ? 1.5 : dUb));
				dLb = -vecdAlpha[i];
				dUb -= vecdAlpha[i];
				dSum += vecdAlpha[i];
			}
			vecdAlpha.back() = 1.0 - dSum;	// 确保alpha的和为1
			
			// 求加权和
			for (auto i = 0; i < m_nGenes; i++)
			{
				double gene = 0.0;
				for (auto j = 0; j < m_nParents; j++)
				{
					gene += vecdAlpha[j] * veciParents[j]->m_vecdChrom[i];
				}
				child.m_vecdChrom.push_back(gene);
			}
			child.m_dFitValue = CalcFitness(child);
			return child;
		}

		void CEvolutionMP::Cross()
		{
			// 执行两次多父体杂交
			m_veciPopulation.push_back(std::move(this->Cross_()));
			m_veciPopulation.push_back(std::move(this->Cross_()));
		}
	}
}