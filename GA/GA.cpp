#include "GA.h"

namespace GA
{ 
	namespace NaiveGA
	{
		// ���캯��
		CEvolution::CEvolution(
			double(*objfunc)(CIndividual&), int n,
			int np, int nc,
			vector<double>&& lb, vector<double>&& ub,
			double pr, double pc
		) :
			m_fcnObj(objfunc), m_nGenerations(n),
			m_nPopSize(np), m_nGenes(nc),
			m_vecdLB(lb), m_vecdUB(ub),
			m_dPr(pr), m_dPc(pc),
			m_fcnConstraints()
		{}

		//��ʼ����Ⱥ
		void CEvolution::Init()
		{

			this->m_veciPopulation.resize(m_nPopSize);
			for (auto i = 0; i < m_nPopSize; i++)
			{
				m_veciPopulation[i].m_vecdChrom.resize(m_nGenes);
				for (auto j = 0; j < m_nGenes; j++)
				{
					//���ֵ��ΧΪ[lb, ub)
					// XXX: ע�⾫��
					double r = RandRange(m_vecdLB[j], m_vecdUB[j]);
					m_veciPopulation[i].m_vecdChrom[j] = r;
				}
				m_veciPopulation[i].m_dFitValue = CalcFitness(m_veciPopulation[i]);	// ������Ӧֵ
			}
			this->m_nCurrGen = 1;	// ���ؼ�����
			_Sort();	// �Ե�ǰ��Ⱥ�и���������򣬴�һ��ʼ��֤������
		}

		//������һ������
		void CEvolution::NextGeneration()
		{
			// ������Ӧֵ��
			// ���ʱ��Ҫ��ȥ���������Ӧֵ
			// ������֤sumֵһ��Ϊ����������������屻ѡ�еĸ���Ϊ0
			// ע�⣺ֻ��ǰm_nPopSize��������ѡ����ֻ��������
			std::for_each(m_veciPopulation.cbegin(),
				m_veciPopulation.cbegin() + m_nPopSize,
				[this](auto e) {m_dSumF += (e.m_dFitValue - m_veciPopulation.back().m_dFitValue); }
			);
			// �Ŵ��������ص�
			for (auto k = 0; k < m_nPopSize; k += 2)
			{
				// ÿ�ε�������2������
				double r = Rand01();
				if (r < m_dPr)
				{
					// ��ֳ
					Reproduce();
				}
				else if (r < (m_dPc + m_dPr))
				{
					// �ӽ�
					Cross();
				}
				else
				{
					// ����
					// ִ�����Σ�ȷ�����������¸���
					Mutate();
					Mutate();
				}
			}

			// ��Ⱥ�ص�
			this->_Sort();	// Ϊ����������
			// ���ڴ�ʱ��Ⱥ��������ֱ���ų�����ĸ���
			m_veciPopulation.erase(m_veciPopulation.begin() + m_nPopSize, m_veciPopulation.end());

			//// ��Ⱥ���ص�
			//// ����ȫ������������Ϊ�Ӵ�
			//m_veciPopulation.erase(m_veciPopulation.begin(), m_veciPopulation.begin()+m_nPopSize);
			//this->_Sort();	// Ϊ��һ������

			++m_nCurrGen;	// ���µ�ǰ����
		}

		// ��ֳ��ѡ�񣩲���
		void CEvolution::Reproduce()
		{
			// ÿ�ζ�����������в���
			const CIndividual& ind1 = Select(), ind2 = Select();
			m_veciPopulation.push_back(ind1);	// ������������ӵ�����β��
			m_veciPopulation.push_back(ind2);
		}

		//�ӽ�����
		void CEvolution::Cross()
		{
			// ���ѡȡ����������Ϊ�����ӽ�
			const CIndividual& father = Select(), mother = Select();//  m_veciPopulation[rand() % m_nPopSize];
			CIndividual child1, child2;

			// �Ը�50%�ĸ��ʽ�����ɢ�ӽ��������ӽ�
			if (Rand01() > 0.5)
			{
				// ������ɢ�ӽ�
				int r = rand() % m_nGenes;
				child1.m_vecdChrom.resize(m_nGenes);
				child2.m_vecdChrom.resize(m_nGenes);
				for (auto i = 0; i < m_nGenes; i++)
				{
					if (i <= r)
					{
						// ����r��ǰ����rλ�ã��Ļ�������
						child1.m_vecdChrom[i] = father.m_vecdChrom[i];
						child2.m_vecdChrom[i] = mother.m_vecdChrom[i];
					}
					else
					{
						// ����r�Ժ�Ļ�������
						child1.m_vecdChrom[i] = father.m_vecdChrom[i];
						child2.m_vecdChrom[i] = mother.m_vecdChrom[i];
					}
				}
			}
			else
			{
				// ���������ӽ�
				double alpha = 0.0;
				for (
					auto itr1 = father.m_vecdChrom.cbegin(), itr2 = mother.m_vecdChrom.cbegin();
					itr1 != father.m_vecdChrom.cend(), itr2 != mother.m_vecdChrom.cend();
					++itr1, ++itr2
					)
				{
					child1.m_vecdChrom.push_back(alpha * (*itr1) + (1.0 - alpha) * (*itr2));
					child2.m_vecdChrom.push_back(alpha * (*itr2) + (1.0 - alpha) * (*itr1));
				}
			}

			child1.m_dFitValue = CalcFitness(child1);
			m_veciPopulation.push_back(std::move(child1));
			child2.m_dFitValue = CalcFitness(child2);
			m_veciPopulation.push_back(std::move(child2));
		}

		// ��������
		void CEvolution::Mutate()
		{
			// ���ѡȡһ��������б���
			// �����и�����Ȼ��ᣬ��һ���������̶ķ�ʽѡ
			//const CIndividual& parent = Select();
			const CIndividual& parent = m_veciPopulation[rand() % m_nPopSize];
			CIndividual child(parent);

			// ÿ��ֻ�仯һ������
			int r = rand() % m_nGenes;
			child.m_vecdChrom[r] = RandRange(m_vecdLB[r], m_vecdUB[r]);

			child.m_dFitValue = CalcFitness(child);
			m_veciPopulation.push_back(std::move(child));
		}

		// ѡ������
		const CIndividual& CEvolution::Select() const
		{
			double minval = m_veciPopulation.back().m_dFitValue;
			if (fabs(m_dSumF) < 1e-15)	return m_veciPopulation[rand() % m_nPopSize];	// ��sumΪ0���������һ������
			double p = Rand01();
			for (auto i = 0; i < m_nPopSize; i++)
			{
				// ���̶�
				p -= (m_veciPopulation[i].m_dFitValue - minval) / m_dSumF;
				if (p < 0)
				{
					return m_veciPopulation[i];
				}
			}
			// ���ǵ��������˴����ؾ���Ⱥ�����һ������
			return m_veciPopulation[m_nPopSize - 1];
		}
	}

	double Rand01()
	{
		return static_cast<double>(rand()) / RAND_MAX;	// ���߾���
	}

	double RandRange(double lb, double ub)
	{
		return Rand01() * (ub - lb) + lb;
	}
}