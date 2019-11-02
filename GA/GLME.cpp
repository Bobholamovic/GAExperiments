#include "GLME.h"
#include <thread>
#include <iostream>
#include <iomanip>
#include <ctime>

using std::endl;

#define VERBOSE
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
			// ʵ���Է����������ݻ�
			// ��һ����ִ���������ݻ�����
			// TODO: ��־�߼��о�Ӧ�ð����ȥ

			clock_t t1 = clock();

			// ȫ���ݻ�
			COUT << "��ʼ����ϣ���ʼȫ���ݻ�" << endl;
			for (auto g = 0; g < m_nGlobalGen__; g++)
			{
				GAMP::CEvolutionMP::NextGeneration();
				COUT << "��" << std::setw(4) << m_nCurrGen << "��\t";
				COUT << std::fixed;
				COUT << "���Ÿ�����Ӧ�ȣ�" << m_veciPopulation[0].m_dFitValue << endl;
			}
			COUT << "ȫ���ݻ���ɣ������� " << m_nGlobalGen__ << " ��" << endl << endl;
			COUT << "�ý׶κ�ʱ: " << std::setprecision(3) << static_cast<double>(clock() - t1) / CLOCKS_PER_SEC << 's' << endl;
			SYSTEM("pause");
			t1 = clock();

			// ���������ӿռ�
			COUT << "���ڹ����ֲ������ӿռ�..." << endl;
			FormClusters();
			int nClusters = m_veciClusters__.size();
			COUT << "������ϣ������ӿռ���ĿΪ " << nClusters << endl << endl;
			COUT << "�ý׶κ�ʱ: " << std::setprecision(3) << static_cast<double>(clock() - t1) / CLOCKS_PER_SEC << 's' << endl;
			SYSTEM("pause");
			t1 = clock();

			// �ֲ��ݻ�
			COUT << "׼��ִ�оֲ��ݻ����������� " << n_threads << " ���߳�" << endl;
			int nLocalGen = m_nGenerations - m_nGlobalGen__;
			int nClustersPerThread = nClusters / n_threads;

			std::unique_ptr<std::thread[]> giThreads(new std::thread[n_threads]);
			// �������̣߳��̼߳�������generationͬ��
			for (auto i = 0; i < n_threads; i++)
			{
				auto itr_s = m_veciClusters__.begin() + std::min(nClustersPerThread * i, nClusters);
				auto itr_t = m_veciClusters__.begin() + std::min(nClustersPerThread * (i+1), nClusters);
				giThreads[i] = std::thread(&CGLME::LocEvolve, itr_s, itr_t, nLocalGen);
			}

			COUT << "����ִ�оֲ��ݻ�����..." << endl;

			// �ȴ������߳��������
			for (auto i = 0; i < n_threads; i++)
			{
				if (giThreads[i].joinable())
					giThreads[i].join();
			}

			COUT << "�ֲ��ݻ���ɣ������� " << nLocalGen << " ��" << endl;
			COUT << "�ý׶κ�ʱ: " << std::setprecision(3) << static_cast<double>(clock() - t1) / CLOCKS_PER_SEC << 's' << endl;
			COUT << endl;

			// ��ʱ���ڸ��ӿռ��ͬ��������ȫ�ֵ�m_nCurrGen�����Ѿ�������
		}

		void CGLME::NextGeneration()
		{
			// ���а汾
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
				// ����Ӧ����Ѱ�Ҵ�����
				// ������K-means����������̣�
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
			// ��ʱ�õ������ĵ�����Ӧ���Ѿ��ǰ�����Ӧֵ�����
			for (auto i = 0; i < nCenters; i++)
			{
				// ���������ӿռ䣨���ɴأ�
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
			veciRets.push_back(pts[0]);	// ��֤��Ӧֵ���Ÿ���һ����ѡȡ
			// XXX: �˴�Ӧ���㷨�Ż�
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
			// �����㷨���صĸ���Ӧ�ñ���ԭ����˳��
			return veciRets;
		}

		vector<CIndividual> CGLME::PickSolutions(double eps)
		{
			// ѡ����ø���
			vector<CIndividual> veciBests;
			for (auto c : m_veciClusters__)
			{
				veciBests.push_back(c.m_veciPopulation[0]);
			}
			// ȥ���ظ����壬�ȼ���Ѱ�Ҿ�������
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