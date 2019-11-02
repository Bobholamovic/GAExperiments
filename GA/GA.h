#pragma once
#include <string>
#include <algorithm>
#include <vector>
#include <functional>
#include <cmath>

using std::vector;

namespace GA
{
	typedef vector<double> CHROM;
	const double MIN_VAL = -1e100;	// 需要根据实际问题调整

	double Rand01();

	double RandRange(double lb, double ub);

	//遗传个体类
	class CIndividual
	{
	public:
		//构造函数
		CIndividual() : m_dFitValue(0.0), m_vecdChrom() {};
		//染色体
		CHROM m_vecdChrom;
		//适应值
		double m_dFitValue;
	};

	using OBJ_FUNC_TYPE=std::function<double(CIndividual&)>;

	namespace NaiveGA
	{
		//进化处理类
		class CEvolution
		{
		protected:
			double m_dSumF = 0.0;	// 当前代所有个体适应值之和，用于加速
		protected:
			static bool Criterion_(const CIndividual& ind1, const CIndividual& ind2)
			{
				return ind1.m_dFitValue > ind2.m_dFitValue;	// 降序排列
			}
			inline void Sort_()
			{
				std::sort(m_veciPopulation.begin(), m_veciPopulation.end(), Criterion_);
			}

		public:
			vector <CIndividual> m_veciPopulation;	// 种群
			int m_nCurrGen;	// 当前代数
			int m_nGenerations;	// 总代数
			const OBJ_FUNC_TYPE m_fcnObj;	// 目标函数
			int m_nPopSize;	// 种群大小
			int m_nGenes;	// 每个染色体上基因数
			double m_dPr;	// 繁殖概率
			double m_dPc;	// 杂交概率
			vector<double> m_vecdLB;	// 变元取值下界
			vector<double> m_vecdUB;	// 变元取值上界
			vector<std::function<double(CHROM&)>>  m_vecfcnConstraints;	// 约束函数

		public:
			//构造函数，调用初始化函数
			CEvolution(
				const OBJ_FUNC_TYPE&, int, int, int, 
				vector<double>&&, vector<double>&&, double, double
			);
			//初始化种群
			void Init();
			//生成下一代
			void NextGeneration();
			// 遗传操作
			// 繁殖算子
			virtual void Reproduce();
			// 交叉算子
			virtual void Cross();
			// 变异算子
			virtual void Mutate();
			// 选择父体
			virtual const CIndividual& Select() const;
			// 适应值函数
			virtual double CalcFitness(CIndividual& ind)
			{
				double dTarVal = this->m_fcnObj(ind);
				for (auto i = 0; i < m_nGenes; i++)
				{
					// 检查每一个分量的取值范围
					if ((ind.m_vecdChrom[i] - m_vecdLB[i]) * (m_vecdUB[i] - ind.m_vecdChrom[i]) < 0)
						return MIN_VAL;
				}

				for (auto fun : this->m_vecfcnConstraints)
				{
					if (fun(ind.m_vecdChrom) < 0)
						return MIN_VAL;
				}
				return dTarVal;
			}
		};
	}

}