#include "GAMP.h"

namespace GA
{
	namespace GAMP
	{
		void CEvolutionMP::_SelectParents(vector<CIndividual const*>& parents)
		{
			// �ȸ������ѡȡ m_nParents ������
			for (auto i = 0; i < m_nParents; i++)
			{
				parents.push_back(&m_veciPopulation[rand()%m_nPopSize]);
			}
		}

		CIndividual CEvolutionMP::_Cross()
		{
			// �ุ���ӽ�
			vector<double> vecdAlpha(m_nParents,0.0);
			CIndividual child;
			double dLb = -0.5, dUb = 1.5;	// ѡȡalpha�����½�
			double dSum = 0.0;
			vector<CIndividual const *> veciParents;

			_SelectParents(veciParents);	// ��ѡ����

			// ������� alpha
			for (auto i = 0; i < m_nParents-1; i++)
			{
				vecdAlpha[i] = RandRange((dLb < -0.5 ? -0.5 : dLb), (dUb > 1.5 ? 1.5 : dUb));
				dLb = -vecdAlpha[i];
				dUb -= vecdAlpha[i];
				dSum += vecdAlpha[i];
			}
			vecdAlpha.back() = 1.0 - dSum;	// ȷ��alpha�ĺ�Ϊ1
			
			// ���Ȩ��
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
			// ִ�����ζุ���ӽ�
			m_veciPopulation.push_back(std::move(this->_Cross()));
			m_veciPopulation.push_back(std::move(this->_Cross()));
		}
	}
}