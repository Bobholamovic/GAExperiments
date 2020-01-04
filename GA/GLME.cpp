#include "GLME.h"
#include <thread>
#include <iostream>
#include <iomanip>
#include <ctime>

using std::endl;

//#define VERBOSE
#ifdef VERBOSE
#define COUT std::cout
#define SYSTEM(x) std::system(x)
#else
#define COUT /##/
#define SYSTEM(x) 
#endif

namespace GA
{
	namespace GLME
	{
		void CGLME::LocEvolve(CGLME::MPE_ITR itr_s, CGLME::MPE_ITR itr_t, int gen)
		{
			for (auto itr = itr_s; itr != itr_t; itr++)
			{
				for (auto g = 0; g < gen; g++)
				{
					itr->NextGeneration();
				}
			}
		}

		void CGLME::ParEvolve(int n_threads)
		{
			// 实验性方案，并行演化
			// 将一次性执行完整个演化过程
			// TODO: 日志逻辑感觉应该剥离出去

			clock_t t1 = clock();

			// 全局演化
			COUT << "初始化完毕，开始全局演化" << endl;
			for (auto g = 0; g < m_nGlobalGen__; g++)
			{
				GAMP::CEvolutionMP::NextGeneration();
				COUT << "第" << std::setw(4) << m_nCurrGen << "代\t";
				COUT << std::fixed;
				COUT << "最优个体适应度：" << m_veciPopulation[0].m_dFitValue << endl;
			}
			COUT << "全局演化完成，共历经 " << m_nGlobalGen__ << " 代" << endl << endl;
			COUT << "该阶段耗时: " << std::setprecision(3) << static_cast<double>(clock() - t1) / CLOCKS_PER_SEC << 's' << endl;
			SYSTEM("pause");
			t1 = clock();

			// 构建搜索子空间
			COUT << "正在构建局部搜索子空间..." << endl;
			FormClusters();
			int nClusters = m_veciClusters__.size();
			COUT << "构建完毕，搜索子空间数目为 " << nClusters << endl << endl;
			COUT << "该阶段耗时: " << std::setprecision(3) << static_cast<double>(clock() - t1) / CLOCKS_PER_SEC << 's' << endl;
			SYSTEM("pause");
			t1 = clock();

			// 局部演化
			COUT << "准备执行局部演化操作，启用 " << n_threads << " 个线程" << endl;
			int nLocalGen = m_nGenerations - m_nGlobalGen__;
			int nClustersPerThread = nClusters / n_threads;

			std::unique_ptr<std::thread[]> giThreads(new std::thread[n_threads]);
			// 启动多线程，线程间无需做generation同步
			for (auto i = 0; i < n_threads; i++)
			{
				auto itr_s = m_veciClusters__.begin() + std::min(nClustersPerThread * i, nClusters);
				auto itr_t = m_veciClusters__.begin() + std::min(nClustersPerThread * (i+1), nClusters);
				giThreads[i] = std::thread(&CGLME::LocEvolve, itr_s, itr_t, nLocalGen);
			}

			COUT << "正在执行局部演化操作..." << endl;

			// 等待所有线程搜索完成
			for (auto i = 0; i < n_threads; i++)
			{
				if (giThreads[i].joinable())
					giThreads[i].join();
			}

			COUT << "局部演化完成，共历经 " << nLocalGen << " 代" << endl;
			COUT << "该阶段耗时: " << std::setprecision(3) << static_cast<double>(clock() - t1) / CLOCKS_PER_SEC << 's' << endl;
			COUT << endl;

			// 此时由于各子空间非同步搜索，全局的m_nCurrGen参数已经无意义
		}

		void CGLME::NextGeneration()
		{
			// 串行版本
			if (m_nCurrGen < m_nGlobalGen__)
			{
				GAMP::CEvolutionMP::NextGeneration();
			}
			else
			{
				if (m_nCurrGen == m_nGlobalGen__)
					this->FormClusters();
				for (auto &c : m_veciClusters__)
				{
					c.NextGeneration();
				}
				m_nCurrGen++;
			}
		}

		void CGLME::FormClusters()
		{
			double epsilon = m_dEpsilon;
			vector<CIndividual> veciCenters;
			int nCenters = 0;

			while (true)
			{
				// 自适应迭代寻找簇中心
				// 或许用K-means代替这个过程？
				veciCenters = GetAllCenters__(m_veciPopulation, epsilon);
				nCenters = veciCenters.size();
				if (nCenters < m_nCentersMin__)
				{
					epsilon /= (1.0 + m_dAlpha);
				}
				else if (nCenters > m_nCentersMax__)
				{
					epsilon *= (1.0 + m_dAlpha);
				}
				else
				{
					break;
				}
			}
			// 此时得到的中心点序列应该已经是按照适应值排序的
			for (auto i = 0; i < nCenters; i++)
			{
				// 定义搜索子空间（生成簇）
				double D = epsilon * ( 1.0 + double(i)/(double(nCenters)-1.0)) / 2.0;
				const auto& center = veciCenters[i];
				CHROM lb(m_nGenes);
				CHROM ub(m_nGenes);
				for (auto j = 0; j < m_nGenes; j++)
				{
					lb[j] = std::max(m_vecdLB[j], center.m_vecdChrom[j] - D);
					ub[j] = std::min(m_vecdUB[j], center.m_vecdChrom[j] + D);
				}
				auto tmpObj = GAMPElite::CEvolutionMPE(
					m_fcnObj, m_nGenerations - m_nGlobalGen__, m_nSubPopSize__,
					m_nGenes, std::move(lb), 
					std::move(ub), m_dPr, m_dPc, m_nParents, 
					m_nElites, m_nCandidates
				);
				tmpObj.m_vecfcnConstraints.insert(tmpObj.m_vecfcnConstraints.end(), 
					m_vecfcnConstraints.begin(), m_vecfcnConstraints.end()
				);
				tmpObj.m_vecfcnConstraints.push_back(Constraint__(D, center.m_vecdChrom));
				tmpObj.Init();
				m_veciClusters__.push_back(std::move(tmpObj));
			}
		}

		vector<CIndividual> CGLME::GetAllCenters__(const vector<CIndividual>& pts, double eps)
		{
			vector<CIndividual> veciRets;
			// veciRets.push_back(pts[rand()%pts.size()]);
			veciRets.push_back(pts[0]);	// 保证适应值最优个体一定被选取
			// XXX: 此处应做算法优化
			bool bFlag = true;
			for (auto itr = pts.begin() + 1; itr != pts.end(); itr++)
			{
				bFlag = true;
				std::for_each(veciRets.begin(), veciRets.end(),
					[=, &bFlag](auto e) {bFlag &= (CalcDist__(e.m_vecdChrom, itr->m_vecdChrom) > eps); }
				);
				if (bFlag)
					veciRets.push_back(*itr);
			}
			// 这种算法返回的个体应该保持原来的顺序
			return veciRets;
		}

		vector<CIndividual> CGLME::PickSolutions(double eps)
		{
			// 选择最好个体
			vector<CIndividual> veciBests;
			for (auto c : m_veciClusters__)
			{
				veciBests.push_back(c.m_veciPopulation[0]);
			}
			// 去除重复个体，等价于寻找聚类中心
			return std::move(GetAllCenters__(veciBests, eps));
		}

		double CGLME::CalcDist__(const CHROM& a, const CHROM& b)
		{
			double sum = 0.0;
			for (unsigned int i = 0; i < a.size(); i++)
				sum += (a[i] - b[i]) * (a[i] - b[i]);
			return sqrt(sum);
		}
	}
}