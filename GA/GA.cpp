#include "GA.h"

namespace GA
{ 
	namespace NaiveGA
	{
		// 构造函数
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

		//初始化种群
		void CEvolution::Init()
		{

			this->m_veciPopulation.resize(m_nPopSize);
			for (auto i = 0; i < m_nPopSize; i++)
			{
				m_veciPopulation[i].m_vecdChrom.resize(m_nGenes);
				for (auto j = 0; j < m_nGenes; j++)
				{
					//随机值范围为[lb, ub)
					// XXX: 注意精度
					double r = RandRange(m_vecdLB[j], m_vecdUB[j]);
					m_veciPopulation[i].m_vecdChrom[j] = r;
				}
				m_veciPopulation[i].m_dFitValue = CalcFitness(m_veciPopulation[i]);	// 更新适应值
			}
			this->m_nCurrGen = 1;	// 拨回计数器
			_Sort();	// 对当前种群中个体进行排序，从一开始保证有序性
		}

		//生成新一代个体
		void CEvolution::NextGeneration()
		{
			// 更新适应值和
			// 求和时需要减去最差个体的适应值
			// 这样保证sum值一定为正，但坏处是最坏个体被选中的概率为0
			// 注意：只从前m_nPopSize个个体中选，即只遍历父代
			std::for_each(m_veciPopulation.cbegin(),
				m_veciPopulation.cbegin() + m_nPopSize,
				[this](auto e) {m_dSumF += (e.m_dFitValue - m_veciPopulation.back().m_dFitValue); }
			);
			// 遗传操作非重叠
			for (auto k = 0; k < m_nPopSize; k += 2)
			{
				// 每次迭代新增2个个体
				double r = Rand01();
				if (r < m_dPr)
				{
					// 繁殖
					Reproduce();
				}
				else if (r < (m_dPc + m_dPr))
				{
					// 杂交
					Cross();
				}
				else
				{
					// 变异
					// 执行两次，确保产生两个新个体
					Mutate();
					Mutate();
				}
			}

			// 种群重叠
			this->_Sort();	// 为新序列排序
			// 由于此时种群向量有序，直接排除靠后的个体
			m_veciPopulation.erase(m_veciPopulation.begin() + m_nPopSize, m_veciPopulation.end());

			//// 种群非重叠
			//// 抛弃全部父代，更换为子代
			//m_veciPopulation.erase(m_veciPopulation.begin(), m_veciPopulation.begin()+m_nPopSize);
			//this->_Sort();	// 为新一代排序

			++m_nCurrGen;	// 更新当前代数
		}

		// 繁殖（选择）操作
		void CEvolution::Reproduce()
		{
			// 每次对两个父体进行操作
			const CIndividual& ind1 = Select(), ind2 = Select();
			m_veciPopulation.push_back(ind1);	// 拷贝副本，添加到向量尾部
			m_veciPopulation.push_back(ind2);
		}

		//杂交算子
		void CEvolution::Cross()
		{
			// 随机选取两个个体作为父本杂交
			const CIndividual& father = Select(), mother = Select();//  m_veciPopulation[rand() % m_nPopSize];
			CIndividual child1, child2;

			// 以各50%的概率进行离散杂交或算术杂交
			if (Rand01() > 0.5)
			{
				// 部分离散杂交
				int r = rand() % m_nGenes;
				child1.m_vecdChrom.resize(m_nGenes);
				child2.m_vecdChrom.resize(m_nGenes);
				for (auto i = 0; i < m_nGenes; i++)
				{
					if (i <= r)
					{
						// 保持r以前（含r位置）的基因序列
						child1.m_vecdChrom[i] = father.m_vecdChrom[i];
						child2.m_vecdChrom[i] = mother.m_vecdChrom[i];
					}
					else
					{
						// 交换r以后的基因序列
						child1.m_vecdChrom[i] = father.m_vecdChrom[i];
						child2.m_vecdChrom[i] = mother.m_vecdChrom[i];
					}
				}
			}
			else
			{
				// 整体算术杂交
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

		// 变异算子
		void CEvolution::Mutate()
		{
			// 随机选取一个父体进行变异
			// 给所有个体均等机会，不一定按照轮盘赌方式选
			//const CIndividual& parent = Select();
			const CIndividual& parent = m_veciPopulation[rand() % m_nPopSize];
			CIndividual child(parent);

			// 每次只变化一个分量
			int r = rand() % m_nGenes;
			child.m_vecdChrom[r] = RandRange(m_vecdLB[r], m_vecdUB[r]);

			child.m_dFitValue = CalcFitness(child);
			m_veciPopulation.push_back(std::move(child));
		}

		// 选择算子
		const CIndividual& CEvolution::Select() const
		{
			double minval = m_veciPopulation.back().m_dFitValue;
			if (fabs(m_dSumF) < 1e-15)	return m_veciPopulation[rand() % m_nPopSize];	// 若sum为0，随机返回一个个体
			double p = Rand01();
			for (auto i = 0; i < m_nPopSize; i++)
			{
				// 轮盘赌
				p -= (m_veciPopulation[i].m_dFitValue - minval) / m_dSumF;
				if (p < 0)
				{
					return m_veciPopulation[i];
				}
			}
			// 考虑到浮点误差，此处返回旧种群中最后一个父体
			return m_veciPopulation[m_nPopSize - 1];
		}
	}

	double Rand01()
	{
		return static_cast<double>(rand()) / RAND_MAX;	// 更高精度
	}

	double RandRange(double lb, double ub)
	{
		return Rand01() * (ub - lb) + lb;
	}
}