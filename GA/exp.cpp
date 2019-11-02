#include <iostream>
#include <iomanip>
#include <ctime>
#include <functional>

#include "GA.h"
#include "GAMP.h"
#include "GAMPElite.h"
#include "GLME.h"

// TODO: 关于输出控制，合并GLME.cpp和此处的宏开关？
#define VERBOSE

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

double Exercise1Proto(GA::CIndividual& ind, int n)
{
	double prod = 1.0;
	const GA::CHROM& x = ind.m_vecdChrom;

	for (auto i = 0; i < n; i++)
	{
		double sum = 0.0;
		for (auto j = 1; j <= 5; j++)
		{
			sum += double(j) * cos((double(j) + 1.0) * x[i] + double(j));
		}
		prod *= sum;
	}

	return -prod;
}

//double Exercise2(GA::CIndividual& ind)
//{
//	double x = ind.m_vecdChrom[0], y = ind.m_vecdChrom[1], z = ind.m_vecdChrom[2];
//
//	double eq1 = sin(x) * cos(y) + sin(y) * cos(z) + sin(z) * cos(x);
//	double eq2 = tan(x) * tan(y) + tan(y) * tan(z) + tan(z) * tan(x);
//	double eq3 = sin(x) * tan(y) * tan(z);
//
//	double v1 = 1.0 - sqrt(6.0) / 2.0;
//	double v2 = 1.0 + 4.0 * sqrt(3.0) / 3.0;
//	double v3 = 1.0 / 2.0;
//
//	return (eq1 - v1) * (v1 - eq1) + (eq2 - v2) * (v2 - eq2) + (eq3 - v3) * (v3 - eq3);
//}

//double Exercise3(GA::CIndividual& ind)
//{
//	// 0-3表示x，4-7表示y
//	auto px = ind.m_vecdChrom.begin(), py = ind.m_vecdChrom.begin()+4;
//
//	double eq1 = *px + *(px + 1) + *(px + 2) + *(px + 3);
//	double eq2 = *px * cos(*py) + *(px + 1) * cos(*(py + 1)) + *(px + 2) * cos(*(py + 2)) + *(px + 3) * cos(*(py + 3));
//	double eq3 = *px * sin(*py) + *(px + 1) * sin(*(py + 1)) + *(px + 2) * sin(*(py + 2)) + *(px + 3) * sin(*(py + 3));
//	double eq4 = *px * cos(*py) * cos(*py) + *(px + 1) * cos(*(py + 1)) * cos(*(py + 1)) + *(px + 2) * cos(*(py + 2)) * cos(*(py + 2)) + \
//		* (px + 3) * cos(*(py + 3)) * cos(*(py + 3));
//	double eq5 = *px * cos(*py) * sin(*py) + *(px + 1) * cos(*(py + 1)) * sin(*(py + 1)) + *(px + 2) * cos(*(py + 2)) * sin(*(py + 2)) + \
//		* (px + 3) * cos(*(py + 3)) * sin(*(py + 3));
//	double eq6 = *px * *px + *(px + 1) * *(px + 1) + *(px + 2) * *(px + 2) + *(px + 3) * *(px + 3);
//	double eq7 = *px * *px * cos(*py) + *(px + 1) * *(px + 1) * cos(*(py+1)) + *(px + 2) * *(px + 2) * cos(*(py + 2)) + \
//		*(px + 3) * *(px + 3) * cos(*(py + 3));
//	double eq8 = *px * *px * sin(*py) + *(px + 1) * *(px + 1) * sin(*(py + 1)) + *(px + 2) * *(px + 2) * sin(*(py + 2)) + \
//		* (px + 3) * *(px + 3) * sin(*(py + 3));
//
//	double v1 = 18048.0;
//	double v2 = 13657.36315298172;
//	double v3 = 5497.052905295088;
//	double v4 = 14508.29635946082;
//	double v5 = 3157.294805334107;
//	double v6 = 105543794.0;
//	double v7 = 91598751.25016867;
//	double v8 = 33470578.99613227;
//
//	return (eq1 - v1) * (v1 - eq1) + (eq2 - v2) * (v2 - eq2) + (eq3 - v3) * (v3 - eq3) + (eq4 - v4) * (v4 - eq4) + (eq5 - v5) * (v5 - eq5) + \
//		(eq6 - v6) * (v6 - eq6) + (eq7 - v7) * (v7 - eq7) + (eq8 - v8) * (v8 - eq8);
//}

void PrintInfo(const GA::CIndividual& ind, int n)
{
	cout << "(";
	for (auto i = 0; i < n; )
	{
		cout << std::fixed;
		cout << std::setprecision(8) << std::setw(12) << ind.m_vecdChrom[i++];
		if (i != n)
			cout << ", ";
	}
	cout << ")\t适应度：" << ind.m_dFitValue << endl;
}

//int main()
//{
//	const double thre = 1e-9;	// 设置阈值
//
//	GA::GAMPElite::CEvolutionMPE ga(
//		Object, 1000, 100, 5, 
//		vector<double>{78, 33, 27, 27, 27}, 
//		vector<double>{102, 45, 45, 45, 45}, 
//		0.0, 1.0, 10, 5, 2
//	);
//
//	// 添加约束项
//	ga.m_vecfcnConstraints.push_back(Constraint1);
//	ga.m_vecfcnConstraints.push_back(Constraint2);
//	ga.m_vecfcnConstraints.push_back(Constraint3);
//
//	//GA::GAMPElite::CEvolutionMPE ga(
//	//	Exercise2, 1000, 100, 3,
//	//	vector<double>{0, 0, 0},
//	//	vector<double>{10, 10, 10},
//	//	0.0, 1.0, 10, 5, 2
//	//);
//
//	//const double PI = 3.1415926;
//	//GA::GAMPElite::CEvolutionMPE ga(
//	//	Exercise3, 100000, 500, 8,
//	//	vector<double>{-1e4, -1e4, -1e4, -1e4, 0, 0, 0, 0},
//	//	vector<double>{1e4, 1e4, 1e4, 1e4, 2*PI, 2*PI, 2*PI, 2*PI},
//	//	0.0, 1.0, 200, 100, 200
//	//);
//
//	srand(time(0));	// 设置随机数种子
//	clock_t t1 = clock();
//	ga.Init();	// 执行初始化
//
//	while (ga.m_nCurrGen <= ga.m_nGenerations)
//	{
//
//#ifdef VERBOSE
//		cout << "第" << std::setw(4) << ga.m_nCurrGen << "代\t";
//		cout << "当前代最好个体：";
//		PrintInfo(ga.m_veciPopulation[0], ga.m_nGenes);
//#endif // VERBOSE
//
//		ga.NextGeneration();
//
//		// 若当前代最好解和最坏解适应值之差小于阈值，停止迭代
//		if (fabs(ga.m_veciPopulation[0].m_dFitValue - ga.m_veciPopulation[ga.m_nPopSize - 1].m_dFitValue) < thre)
//			break;
//
//	}
//
//	cout << "耗时: " << std::setprecision(3) << static_cast<double>(clock() - t1) / CLOCKS_PER_SEC << 's' << endl;
//	cout << "历经代数：" << ga.m_nCurrGen-1 << endl;
//	cout << "历史最佳个体：";
//	PrintInfo(ga.m_veciPopulation[0], ga.m_nGenes);
//
//	return 0;
//}

int main()
{
	//int nGlobalGen = 1000;

	srand(time(0));	// 设置随机数种子

	//GA::GLME::CGLME ga(
	//	Object, 1000, 300, 5, 
	//	GA::CHROM{78, 33, 27, 27, 27}, 
	//	GA::CHROM{102, 45, 45, 45, 45}, 
	//	0.0, 1.0, 10, 5, 2, nGlobalGen
	//);

	//// 添加约束项
	//ga.m_vecfcnConstraints.push_back(Constraint1);
	//ga.m_vecfcnConstraints.push_back(Constraint2);
	//ga.m_vecfcnConstraints.push_back(Constraint3);

	using std::placeholders::_1;
	int n = 3;
	auto exercise1 = std::bind(Exercise1Proto, _1, n);
	GA::GLME::CGLME ga(
		exercise1,
		3500, 1000, n, 
		GA::CHROM(n, -10), 
		GA::CHROM(n, 10),
		0.0, 1.0, 20, 10, 2, 3000
	);

	// 配置参数
	ga.ConfigLocalParams(80, 160, 0.1, 0.2, 100);

	clock_t t1 = clock();
	ga.Init();	// 执行初始化

//	while (ga.m_nCurrGen <= ga.m_nGenerations)
//	{
//
//#ifdef VERBOSE
//		if (ga.m_nCurrGen < nGlobalGen)
//		{
//			cout << "第" << std::setw(4) << ga.m_nCurrGen << "代\t";
//			cout << "当前代最好个体：";
//			PrintInfo(ga.m_veciPopulation[0], ga.m_nGenes);
//		}
//		else
//		{
//			cout << "第" << std::setw(4) << ga.m_nCurrGen << "代\n";
//		}
//
//#endif // VERBOSE
//
//		ga.NextGeneration();
//
//		//// 若当前代最好解和最坏解适应值之差小于阈值，停止迭代
//		//if (fabs(ga.m_veciPopulation[0].m_dFitValue - ga.m_veciPopulation[ga.m_nPopSize - 1].m_dFitValue) < thre)
//		//	break;
//
//	}

	ga.ParEvolve(16);

	vector<GA::CIndividual> solutions = ga.PickSolutions(1e-3);
	cout << "总过程耗时: " << std::setprecision(3) << static_cast<double>(clock() - t1) / CLOCKS_PER_SEC << 's' << endl << endl;
	cout << "共找到 " << solutions.size() << " 个可能解，分别为：" << endl;
	system("PAUSE");

	// 对解按适应值进行排序
	std::sort(
		solutions.begin(), solutions.end(),
		[](const GA::CIndividual& a, const GA::CIndividual& b) {return a.m_dFitValue > b.m_dFitValue; }
	);
	auto counter = 0;
	for (const auto& s : solutions)
	{
		cout << ++counter << '\t';
		PrintInfo(s, ga.m_nGenes);
	}
		

	return 0;
}