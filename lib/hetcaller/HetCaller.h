// ***************************************************************************
// HetCaller.h (c) 2020 Zhenhua Yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#ifndef _HETCALLER_H
#define _HETCALLER_H

#include <vector>
#include <map>
#include <pthread.h>

#include "Matrix.h"

class Solution {
	public:
		Solution() {}
		Solution(int num_cluster, double alpha, double beta, Matrix<int>& states) :
			num_cluster(num_cluster), alpha(alpha), beta(beta), states(states) {}
		
		int num_cluster; //number of clusters
		bool valid;
		double ll; //log-likelihood
		double bic; //Bayesian information criterion
		double acc; //accuracy
		double alpha; //false positive rate
		double beta; //false negative rate
		Matrix<int> states; //states of clusters
		vector<uint32_t> indices;
};

class HetCaller {
	private:
		Matrix<int> obsData;
		Matrix<int> data_processed;
		Matrix<float> data_for_pca;
		Matrix<int> realData;
		vector<string> cLabels;
		vector<string> mLabels;
		bool homo_muta;
		double missing_rate;
		
		vector<int> doublet_indxs;
		
		//k-means
		vector<vector<float>> data_for_cluster;
		float min_dist; 
		vector<uint32_t> best_indices;
		map<int, int> cell_assignments;
		
		Matrix<int> para_updates;
		vector<Solution> candi_s;
		vector<Solution> candi_s1;
		int best_s_indx;
		
		Matrix<int> clonal_tree;
		
		pthread_mutex_t pm;
		
		void preProcess();
		
		static void* clusterMutations(const void *arg);
		void saveClusterResults(float dist, vector<uint32_t>& indices);
		
		void predict(int num_cluster);
		static void* inferSolution(const void *arg);
		//static void* inferParas(const void *arg);
		void inferParas(Solution& s, const vector<int>& para_updates);
		
		double getObsProb(int gt, int observed, double alpha, double beta);
		void saveSolution(Solution& s);
		
		void printSolution(Solution& s);
		void validity(Solution& s);
		
		void loadMutationData();
		void loadRealData();
		void loadCellLabels();
		void loadMutaLabels();
		
		void evalAccuracy(Solution& s);
		void buildTree();
		void saveResults();
		void saveTree();
		
	public:
		HetCaller();
	
		Matrix<int>& getObsData() {return obsData;}
		Matrix<int>& getProcessedData() {return data_processed;}
		vector<vector<float>>& getDataForCluster() {return data_for_cluster;}
		vector<uint32_t>& getBestIndices() {return best_indices;}
		
		Solution& getSolution(int i) {return candi_s[i];}
		Matrix<int>& getParaUpdates() {return para_updates;}
		
		bool getHomoMuta() {return homo_muta;}
		
		void loadData();
		void call();
};

#endif
