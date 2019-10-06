#include <iostream>
#include <iomanip>
#include <time.h>

#include "GA.h"
#include "GAMP.h"
#include "GAMPElite.h"

using std::cout;
using std::endl;

// 目标函数
double Object(GA::CIndividual& ind)
{
	GA::CHROM& x = ind.m_vecdChrom;
	return -(5.3578547 * x[2] * x[2] + 0.8356891 * x[0] * x[4] + 37.293239 * x[0] - 40792.141);
}

double Constraint1(GA::CHROM& x)
{
	double val = 85.334407 + 0.0056858 * x[1] * x[4] + 0.0006262*x[0]*x[3] - 0.0022053 * x[2] * x[4];
	return val * (92.0 - val);
}

double Constraint2(GA::CHROM& x)
{
	double val = 80.51249 + 0.0071317 * x[1] * x[4] + 0.0029955 * x[0] * x[1] + 0.0021813 * x[2]*x[2];
	return (val-90.0) * (110.0 - val);
}

double Constraint3(GA::CHROM& x)
{
	double val = 9.300961 + 0.0047026 * x[2] * x[4] + 0.0012547 * x[0] * x[2] + 0.0019085 * x[2] * x[3];
	return (val - 20.0) * (25.0 - val);
}

int main()
{
	const double thre = 1e-14;

	GA::GAMPElite::CEvolutionMPE ga(
		Object, 1000, 1000, 5, 
		vector<double>{78, 33, 27, 27, 27}, 
		vector<double>{102, 45, 45, 45, 45}, 
		0.0, 1.0
	);

	// 添加约束项
	ga.m_fcnConstraints.push_back(Constraint1);
	ga.m_fcnConstraints.push_back(Constraint2);
	ga.m_fcnConstraints.push_back(Constraint3);

	srand(time(0));	// 设置随机数种子
	ga.Init();	// 执行初始化

	while (ga.m_nCurrGen <= ga.m_nGenerations)
	{
		cout << "第" << std::setw(4) << ga.m_nCurrGen << "代\t";
		cout << "当前代最好个体: (";
		for (auto i = 0; i < ga.m_nGenes; )
		{
			cout << std::setprecision(8) << std::setw(8) << ga.m_veciPopulation[0].m_vecdChrom[i++];
			if (i != ga.m_nGenes)
				cout << ", ";
		}
		cout << ")\t适应度: " << std::setprecision(8) << std::setw(8) << ga.m_veciPopulation[0].m_dFitValue << endl;

		// 若当前代最好解和最坏解适应值之差小于阈值，停止迭代
		if (fabs(ga.m_veciPopulation[0].m_dFitValue - ga.m_veciPopulation[ga.m_nPopSize - 1].m_dFitValue) < thre)
			break;

		ga.NextGeneration();
	}

	return 0;
}
