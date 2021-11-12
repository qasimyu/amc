// ***************************************************************************
// HetCaller.cpp (c) 2021 zhenhua yu <qasim0208@163.com>
// Health Informatics Lab, Ningxia University
// All rights reserved.

#include <iostream>
#include <cstdio>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <limits>
#include <array>
#include <tuple>

#include "HetCaller.h"
#include "split.h"
#include "MyDefine.h"
#include "PCA.h"
#include "dkm.hpp"
#include "fastcluster.h"

using namespace std;

HetCaller::HetCaller() {
	pthread_mutex_init(&pm, NULL);
}

void HetCaller::loadData() {
	loadMutationData();
	loadRealData();
	loadCellLabels();
	loadMutaLabels();
}

void HetCaller::loadMutationData() {
	string inputFile = config.getStringPara("input");
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	
	FILE *fp;
	char buf[500];
	string cmd = "cat "+inputFile+" | wc -l";
	fp = popen(cmd.c_str(), "r");
	if(!fp) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	char *p = fgets(buf, 500, fp);
	fclose(fp);
	int num_muta = atoi(buf);
	
	int i, j, k, indx;
	long line_num = 0;
	string line;
	int num_cell = 0;
	homo_muta = false;
	while(getline(ifs, line)) {
		line_num++;
		vector<string> fields = split(line, '\t');
		if(num_cell == 0) {
			num_cell = fields.size();
			obsData.resize(num_cell, num_muta, false);
		}
		if(num_cell != fields.size()) {
			cerr << "Error: malformed input file " << inputFile <<
                    ", inconsistent number of fields @line " << line_num << endl;
            cerr << buf << endl;
			exit(1);
		}
		for(i = 0; i < num_cell; i++) {
			indx = i*num_muta+line_num-1;
			k = atoi(fields[i].c_str());
			if(k == 2) {
				homo_muta = true;
			}
			obsData[indx] = k;
		}
	}
	ifs.close();
	assert(line_num == num_muta);
	
	cerr << "Total " << num_cell << " cells and " << num_muta << " mutations were loaded from file " << inputFile << endl;
}

void HetCaller::loadRealData() {
	string inputFile = config.getStringPara("rinput");
	if(inputFile.empty()) {
		return;
	}
	
	FILE *fp;
	char buf[500];
	string cmd = "cat "+inputFile+" | wc -l";
	fp = popen(cmd.c_str(), "r");
	if(!fp) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	char *p = fgets(buf, 500, fp);
	fclose(fp);
	int num_muta = atoi(buf)-1;
	
	int i, j, k, indx;
	long line_num = 0;
	string line;
	int num_cell = 0;
	
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	getline(ifs, line);
	vector<string> fields = split(line, ':');
	if(fields.size() > 1) {
		fields = split(fields[1], '\t');
		for(i = 0; i < fields.size(); i++) {
			doublet_indxs.push_back(atoi(fields[i].c_str())-1);
		}
	}
	
	while(getline(ifs, line)) {
		line_num++;
		vector<string> fields = split(line, '\t');
		if(num_cell == 0) {
			num_cell = fields.size();
			realData.resize(num_cell, num_muta, false);
		}
		if(num_cell != fields.size()) {
			cerr << "Error: malformed input file " << inputFile <<
                    ", inconsistent number of fields @line " << line_num << endl;
            cerr << buf << endl;
			exit(1);
		}
		for(i = 0; i < num_cell; i++) {
			indx = i*num_muta+line_num-1;
			k = atoi(fields[i].c_str());
			realData[indx] = k;
		}
	}
	ifs.close();
	assert(line_num == num_cell);
	
	cerr << "real mutation data were loaded from file " << inputFile << endl;
}

void HetCaller::loadCellLabels() {
	string inputFile = config.getStringPara("clabel");
	if(inputFile.empty()) {
		char buf[100];
		for(int i = 0; i < obsData.getROWS(); i++) {
			sprintf(buf, "%d", i+1);
			cLabels.push_back(buf);
		}
		return;
	}
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	
	long line_num = 0;
	string line;
	while(getline(ifs, line)) {
		line_num++;
		if(line.empty()) {
			cerr << "Error: malformed label file " << inputFile <<
                    ", empty item is not allowed @line " << line_num << endl;
			exit(1);
		}
		cLabels.push_back(line);
	}
	ifs.close();
	
	if(cLabels.size() != obsData.getROWS()) {
		cerr << "Error: the number of cell labels is not consistent with the number of cells." << endl;
		exit(1);
	}
	
	cerr << "the labels of cells were loaded from file " << inputFile << endl;
}

void HetCaller::loadMutaLabels() {
	string inputFile = config.getStringPara("mlabel");
	if(inputFile.empty()) {
		char buf[100];
		for(int i = 1; i <= obsData.getCOLS(); i++) {
			sprintf(buf, "%d", i);
			mLabels.push_back(buf);
		}
		return;
	}
	ifstream ifs;
	ifs.open(inputFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: cannot open file " << inputFile << endl;
		exit(-1);
	}
	
	long line_num = 0;
	string line;
	while(getline(ifs, line)) {
		line_num++;
		if(line.empty()) {
			cerr << "Error: malformed label file " << inputFile <<
                    ", empty item is not allowed @line " << line_num << endl;
			exit(1);
		}
		mLabels.push_back(line);
	}
	ifs.close();
	
	if(mLabels.size() != obsData.getCOLS()) {
		cerr << "Error: the number of mutation labels is not consistent with the number of mutations." << endl;
		exit(1);
	}
	
	cerr << "the labels of mutations were loaded from file " << inputFile << endl;
}

void HetCaller::preProcess() {
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	int i, j, k;
	
	data_processed.resize(num_muta, num_cell, false);
	data_for_pca.resize(num_muta, num_cell, false);
	missing_rate = 0;
	for(i = 0; i < num_cell; i++) {
		for(j = 0; j < num_muta; j++) {
			k = obsData[i*num_muta+j];
			if(homo_muta) {
				if(k == 3) {
					data_for_pca[j*num_cell+i] = 0.5;
					data_processed[j*num_cell+i] = 0;
					missing_rate++;
				}
				else if(k == 2) {
					data_for_pca[j*num_cell+i] = 0.99;
					data_processed[j*num_cell+i] = 1;
				}
				else {
					data_for_pca[j*num_cell+i] = k;
					data_processed[j*num_cell+i] = k;
				}
			}
			else {
				if(k == 3) {
					data_for_pca[j*num_cell+i] = 0.5;
					data_processed[j*num_cell+i] = 0;
					missing_rate++;
				}
				else {
					data_for_pca[j*num_cell+i] = k;
					data_processed[j*num_cell+i] = k;
				}
			}
		}
	}
	missing_rate /= (num_cell*num_muta);
}

void HetCaller::call() {
	int i, j, k, n;
	int num_muta = obsData.getCOLS();
	int num_cell = obsData.getROWS();
	
	//get data for PCA analysis
	preProcess();
	
	Matrix<float> reduced_data;
	//if(num_muta >= num_cell) {
		//perform PCA, keeping 90% energy
		PCA pca(data_for_pca);
		//pca.first_k_ONB(num_muta-1);
		pca.important_ONB(0.9);
		reduced_data = pca.getProj();
	//}
	//else {
		//reduced_data = data_for_pca;
	//}
	/*
	for(i = 0; i < 2; i++) {
		reduced_data.Row(i).Print();
		cerr << endl;
	}
	*/
	int num_v = reduced_data.getCOLS();
	cerr << "Dimension reduces to " << reduced_data.getROWS() << "X" << reduced_data.getCOLS() << endl;
		
	/*** group mutations into clusters ***/
	
	
	for(i = 0; i < reduced_data.getROWS(); i++) {
		vector<float> tmp;
		tmp.reserve(reduced_data.getCOLS());
		for(j = 0; j < reduced_data.getCOLS(); j++) {
			tmp.push_back(reduced_data[i*reduced_data.getCOLS()+j]);
		}
		data_for_cluster.push_back(tmp);
	}
	
	
	/*** find the optimal number of clusters based on BIC ***/
	int maxc = config.getIntPara("maxc");
	if(maxc == -1) {
		maxc = min(100, num_muta);
	}
	else {
		maxc = min(maxc, num_muta);
	}
	int km_reps = config.getIntPara("km_reps");
	
	double min_bic = numeric_limits<double>::max();
	int count = 0;
	for(int num_cluster = 1; num_cluster <= maxc; num_cluster++) {
		/*** group mutations into clusters ***/
		min_dist = numeric_limits<float>::max();
		vector<int*> t_paras;
		for(i = 0; i < km_reps; i++) {
			int* paras = new int[2];
			paras[0] = num_cluster; paras[1] = i;
			t_paras.push_back(paras);
			threadpool->pool_add_work(&HetCaller::clusterMutations, paras, i);
		}
		threadpool->wait();
		for(i = 0; i < km_reps; i++) {
			delete[] t_paras[i];
		}
		
		/*** reasoning the mutation states of each cluster ***/
		predict(num_cluster);
		candi_s[num_cluster-1].indices = best_indices;
		
		cerr << "--------------- screening report -----------------" << endl;
		printSolution(candi_s[num_cluster-1]);
		cerr << "--------------- screening report -----------------" << endl;
		
		
		if(candi_s[num_cluster-1].valid && candi_s[num_cluster-1].bic < min_bic) {
			min_bic = candi_s[num_cluster-1].bic;
			best_s_indx = num_cluster-1;
			count = 0;
		}
		else {
			count++;
			if(count >= 10) {
				break;
			}
		}
	}
	
	best_s_indx = 0;
	for(i = 1; i < candi_s.size(); i++) {
		if(candi_s[i].bic < candi_s[best_s_indx].bic) {
			best_s_indx = i;
		}
	}
	
	best_indices = candi_s[best_s_indx].indices;
	
	/*** best solution ***/
	cerr << "--------------- best solution -----------------" << endl;
	printSolution(candi_s[best_s_indx]);
	cerr << "-----------------------------------------------" << endl;
	
	saveResults();
	
	string scite = config.getStringPara("scite");
	if(scite.empty()) {
		return;
	}
	/*** construct subclonal tree using SCITE***/
	buildTree();
	
	/*** save phylogenetic tree ***/
	saveTree();
	
}

void HetCaller::predict(int num_cluster) {
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	int i, j, epoches;

	if(num_cluster == num_muta) {
		epoches = 1;
	}
	else {
		int N = max(1,num_muta/num_cluster);
		double tmp = 1.0;
		for(i = num_cluster+1, j = 1; i <= num_muta; i++, j++) {
			tmp *= 1.0*i/(i-num_cluster);
			if(j <= num_cluster) {
				tmp /= N;
			}
		}
		while(j <= num_cluster) {
			tmp /= N;
			j++;
		}
		epoches = round(tmp);
		if(epoches > 100) {
			epoches = 100;
		}
		if(epoches < 50) {
			epoches = 50;
		}
	}
	
	candi_s1.clear();
	for(i = 0; i < epoches; i++) {
		threadpool->pool_add_work(&HetCaller::inferSolution, &num_cluster, i);
	}
	threadpool->wait();
	
	int best_indx = -1;
	for(i = 0; i < candi_s1.size(); i++) {
		if(!candi_s1[i].valid) {
			continue;
		}
		if(best_indx == -1) {
			best_indx = i;
		}
		
		else if(candi_s1[i].ll > candi_s1[best_indx].ll) {
			best_indx = i;
		}
	}
	if(best_indx == -1) {
		best_indx = 0;
		for(i = 1; i < candi_s1.size(); i++) {
			if(candi_s1[i].ll > candi_s1[best_indx].ll) {
				best_indx = i;
			}
		}
	}
	
	evalAccuracy(candi_s1[best_indx]);
	double lambda = min(0.15, 50.0/sqrt(num_muta*num_cell));
	candi_s1[best_indx].bic = -2*candi_s1[best_indx].ll+lambda*log(num_muta)*(num_cluster*num_cell+2);
	candi_s1[best_indx].num_cluster = num_cluster;
	candi_s.push_back(candi_s1[best_indx]);
}

void* HetCaller::inferSolution(const void *arg) {
	int num_cluster = *((int*) arg);
	int i, j, k, n;
	
	Matrix<int>& data = hetcaller.getProcessedData();
	vector<uint32_t>& indices = hetcaller.getBestIndices();
	
	int num_cell = data.getCOLS();
	int num_muta = data.getROWS();
	
	Matrix<int> states(num_cluster, num_cell);
	
	map<int, vector<int>> muta_assignments;
	for(i = 0; i < indices.size(); i++) {
		muta_assignments[indices[i]].push_back(i);
	}
	
	for(i = 0; i < num_cluster; i++) {
		vector<int>& muta_indxs = muta_assignments[i];
		k = threadpool->randomInteger(0, muta_indxs.size());
		k = muta_indxs[k];
		for(j = 0; j < num_cell; j++) {
			states[i*num_cell+j] = data[k*num_cell+j];
		}
	}
	
	// search for optimal initial values of beta
	double alpha = config.getRealPara("alpha");
	double beta = config.getRealPara("beta");
	
	vector<double> alphas, betas;
	bool alpha_fixed, beta_fixed;
	if(alpha >= 0) {
		alphas.push_back(alpha);
		alpha_fixed = true;
	}
	else {
		alphas.push_back(0.01);
		alpha_fixed = false;
	}
	if(beta >= 0) {
		betas.push_back(beta);
		beta_fixed = true;
	}
	else {
		betas.push_back(0.01);
		beta_fixed = false;
	}
	
	vector<int> para_updates(3, 1);
	para_updates[0] = !alpha_fixed;
	para_updates[1] = !beta_fixed;
	
	//record best solution
	Solution best_s;
	best_s.ll = numeric_limits<long>::min();
	
	for(i = 0; i < alphas.size(); i++) {
		for(j = 0; j < betas.size(); j++) {
			Solution s(num_cluster, alphas[i], betas[j], states);
			hetcaller.inferParas(s, para_updates);
			if(s.ll > best_s.ll) {
				best_s = s;
			}
		}
	}
	//best_s.indices = indices;
	//hetcaller.calculateInterScore(best_s);
	best_s.valid = true;
	hetcaller.saveSolution(best_s);
	
	return NULL;
}

void HetCaller::inferParas(Solution& s, const vector<int>& para_updates) {
	int i, j, k, n;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	
	int num_cluster = s.num_cluster;
	double& alpha = s.alpha;
	double& beta = s.beta;
	
	Matrix<int>& states = s.states;
	Matrix<int> states_n(num_cluster, num_cell);
	
	double max_alpha = 0.99, max_beta = 0.99;
	double eps = numeric_limits<long double>::epsilon();
	
	double ll, pre_ll = numeric_limits<long>::min();
	int iter = 0, max_iter = config.getIntPara("max_iter");
	while(iter < max_iter) {
		/*** calculate log-likelihood ***/
		ll = 0;
		for(i = 0; i < num_cell; i++) {
			for(j = 0; j < num_muta; j++) {
				k = best_indices[j];
				double prob = getObsProb(states[k*num_cell+i], obsData[i*num_muta+j], alpha, beta);
				ll += log(prob);
			}
		}
		
		// update alpha and beta
		double tmp1 = 0, tmp2 = 0, tmp3 = 0, tmp4 = 0, tmp5 = 0;
		for(i = 0; i < num_cell; i++) {
			for(j = 0; j < num_muta; j++) {
				if(obsData[i*num_muta+j] == 3) {
					continue;
				}
				k = best_indices[j];
				int s = obsData[i*num_muta+j];
				int a = (1-s)*(2-s)/2;
				int b = s*(2-s);
				int c = s*(s-1)/2;
				if(homo_muta) {
					tmp1 += alpha*(a+c+b*states[k*num_cell+i]);
					tmp2 += -2*c-2*states[k*num_cell+i]*a-2*states[k*num_cell+i]*b;
					tmp3 += 2*c*(1-alpha)+2*states[k*num_cell+i]*a*(1-alpha);
					tmp4 += (1-states[k*num_cell+i])*(a+b+c)*(2-beta);
					tmp5 += 2*(1-states[k*num_cell+i])*(b+c);
				}
				else {
					tmp1 += (1-states[k*num_cell+i])*s;
					tmp2 += (1-states[k*num_cell+i]);
					tmp3 += states[k*num_cell+i]*(1-s);
					tmp4 += states[k*num_cell+i];
				}
			}
		}
		double alpha_n, beta_n;
		if(homo_muta) {
			if(alpha <= eps) {
				beta_n = -tmp3/tmp2;
			}
			else {
				beta_n = 0.5*(-tmp2-sqrt(tmp2*tmp2-4*tmp1*tmp3))/(tmp1+eps);
			}
			alpha_n = tmp5/tmp4;
		}
		else {
			alpha_n = tmp1/(tmp2+eps);
			beta_n = tmp3/(tmp4+eps);
		}
		
		// update states
		if(para_updates[2] == 1) {
			for(k = 0; k < num_cluster; k++) {
				for(i = 0; i < num_cell; i++) {
					double f1 = 0, f2 = 0;
					for(j = 0; j < num_muta; j++) {
						if(best_indices[j] != k) {
							continue;
						}
					
						if(obsData[i*num_muta+j] == 3) {
							continue;
						}
						int s = obsData[i*num_muta+j];
						if(homo_muta) {
							f1 += pow(1-s, 2)*log(beta/2+eps)+s*(2-s)*log(1-beta+eps);
							f2 += (1-s)*(2-s)/2*log(1-alpha-alpha*beta/2+eps)+s*(2-s)*log(alpha+eps)+s*(s-1)/2*log(alpha*beta/2+eps);
						}
						else {
							f1 += (s*log(1-beta+eps)+(1-s)*log(beta+eps));
							f2 += (s*log(alpha+eps)+(1-s)*log(1-alpha+eps));
						}
						
					}
					states_n[k*num_cell+i] = (f1 > f2)? 1:0;
				}
			}
			states = states_n;
		}
		
		if(para_updates[0]) {
			alpha = isnan(alpha_n)? alpha:alpha_n;
			alpha = (alpha > max_alpha)? max_alpha:alpha;
		}
		if(para_updates[1]) {
			beta = isnan(beta_n)? beta:beta_n;
			beta = (beta > max_beta)? max_beta:beta;
		}
		
		if(fabs(pre_ll-ll) < 1e-4) {
			break;
		}
		else {
			pre_ll = ll;
		}
		iter++;
	}
	s.ll = ll;
}

void HetCaller::printSolution(Solution& s) {
	cerr << "#clusters: " << s.num_cluster << endl;
	cerr << "LL: " << s.ll << endl;
	cerr << "BIC: " << s.bic << endl;
	//cerr << "valid: " << s.valid << endl;
	
	if(s.acc >= 0) {
		cerr << "Acc: " << s.acc << endl;
	}
	cerr << "alpha: " << s.alpha << endl;
	cerr << "beta: " << s.beta << endl;
}

void HetCaller::validity(Solution& s) {
	int i, j;
	int num_cluster = s.num_cluster; 
	Matrix<int>& states = s.states;
	int num_cell = states.getCOLS();
	s.valid = true;
	for(i = 0; i < num_cluster; i++) {
		for(j = 0; j < num_cell; j++) {
			if(states[i*num_cell+j] == 1) {
				break;
			}
		}
		if(j == num_cell) {
			s.valid = false;
			break;
		}
	}
}

void* HetCaller::clusterMutations(const void *arg) {
	int* tmp = (int*) arg;
	int num_cluster = tmp[0];
	int seed = tmp[1];
	
	vector<vector<float>>& data_for_cluster = hetcaller.getDataForCluster();
	
	float dist, score;
	dkm::clustering_parameters<float> paras(num_cluster);
	//paras.set_random_seed(threadpool->randomInteger(0, 10000));
	paras.set_random_seed(seed);
	paras.set_max_iteration(100);
	auto cluster_data = dkm::kmeans_lloyd(data_for_cluster, paras, dist, score);
	
	vector<uint32_t>& indices = get<1>(cluster_data);
	Matrix<int> flags(1, num_cluster, 0);
	for(int i = 0; i < indices.size(); i++) {
		flags[indices[i]] = 1;
	}
	if(flags.sum() == num_cluster) {
		hetcaller.saveClusterResults(dist, get<1>(cluster_data));
	}
	
	return NULL;
}

void HetCaller::saveClusterResults(float dist, vector<uint32_t>& indices) {
	pthread_mutex_lock(&pm);
	
	if(min_dist > dist) {
		min_dist = dist;
		best_indices = indices;
	}
	pthread_mutex_unlock(&pm);
}

void HetCaller::evalAccuracy(Solution& s) {
	if(realData.getROWS() == 0) {
		s.acc = -1;
		return;
	}
	int i, j, k;
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	
	Matrix<int>& states = s.states;
	double acc = 0;
	int n = 0;
	for(i = 0; i < num_cell; i++) {
		for(j = 0; j < doublet_indxs.size(); j++) {
			if(doublet_indxs[j] == i) {
				break;
			}
		}
		if(j < doublet_indxs.size()) {
			continue;
		}
		n++;
		for(j = 0; j < num_muta; j++) {
			k = best_indices[j];
			k = states[k*num_cell+i];
			if(k == realData[i*num_muta+j])	acc++;
		}
	}
	s.acc = acc/(n*num_muta);
}

double HetCaller::getObsProb(int gt, int observed, double alpha, double beta) {
	double prob = 1, eps = numeric_limits<double>::epsilon();
	if(gt == 1) {
		if(observed == 0) {
			prob = (homo_muta)? beta/2+eps:beta+eps;
		}
		else if(observed == 1) {
			prob = 1-beta+eps;
		}
		else if(observed == 2) {
			prob = beta/2+eps;
		}
	}
	else {
		if(observed == 0) {
			prob = (homo_muta)? 1-alpha-alpha*beta/2+eps:1-alpha+eps;
		}
		else if(observed == 1) {
			prob = alpha+eps;
		}
		else if(observed == 2) {
			prob = alpha*beta/2+eps;
		}
	}
	return prob;
}

void HetCaller::saveSolution(Solution& s) {
	pthread_mutex_lock(&pm);
	candi_s1.push_back(s);
	pthread_mutex_unlock(&pm);
}

void HetCaller::buildTree() {
	string scite = config.getStringPara("scite");
	if(scite.empty()) {
		return;
	}
	
	char buf[5000];
	int num_cluster = candi_s[best_s_indx].num_cluster;
	string outputPrefix = config.getStringPara("output");
	
	//sprintf(buf, "%s_clusters_%d.scite_input", outputPrefix.c_str(), num_cluster);
	//string dataFile = buf;
	string dataFile = outputPrefix+".scite_input";
	
	int num_cell = obsData.getROWS();
	
	sprintf(buf, "%s -i %s -n %d -m %d", scite.c_str(), dataFile.c_str(), num_cluster, num_cell);
	string cmd = buf;
	cmd = cmd+" -r 1 -l 50000 -fd 5e-6 -ad 5e-6 -e 0 -a -max_treelist_size 1";
	sprintf(buf, "%s -o %s > %s.scite.log 2>&1", cmd.c_str(), outputPrefix.c_str(), outputPrefix.c_str());
	cmd = buf;
	
	//cerr << cmd << endl;
	system(cmd.c_str());
	
	
	clonal_tree.resize(1, num_cluster+1, false);
	clonal_tree.set(-1);
	
	string gvFile = outputPrefix+"_ml0.gv";
	ifstream ifs;
	ifs.open(gvFile.c_str());
	if(!ifs.is_open()) {
		cerr << "Error: fail to fetch results of SCITE, file not found: " << gvFile << endl;
		exit(-1);
	}
	
	int i, j, k;
	long line_num = 0;
	string line;
	getline(ifs, line);getline(ifs, line);
	
	while(getline(ifs, line)) {
		i = line.find("node");
		if(i != string::npos) {
			break;
		}
		i = line.find("->");
		string tmp = trim(line.substr(0, i));
		j = atoi(tmp.c_str());
		tmp = trim(line.substr(i+2, line.length()-i-2));
		k = atoi(tmp.c_str());
		//cerr << j << "->" << k << endl;
		clonal_tree[k-1] = j-1;
	}
	
	while(getline(ifs, line)) {
		i = line.find("->");
		if(i == string::npos) {
			break;
		}
		string tmp = trim(line.substr(0, i));
		j = atoi(tmp.c_str());
		tmp = trim(line.substr(i+2, line.length()-i-2));
		k = atoi(tmp.substr(1).c_str());
		//cerr << j << "->" << k << endl;
		cell_assignments[k] = j-1;
	}
	
	ifs.close();
	
	cmd  = "rm "+outputPrefix+".samples";
	system(cmd.c_str());
	
}

void HetCaller::saveResults() {
	int i, j, k;
	ofstream ofs;
	string fn, outputPrefix = config.getStringPara("output");
	
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	/*
	fn = outputPrefix+".solution.txt";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	ofs << "alpha = " << candi_s[best_s_indx].alpha << endl;
	ofs << "beta = " << candi_s[best_s_indx].beta << endl;
	//ofs << "K = " << candi_s[best_s_indx].num_cluster << endl;
	ofs.close();
	*/
	
	int num_cluster = candi_s[best_s_indx].num_cluster;
	Matrix<int>& states = candi_s[best_s_indx].states;
	
	fn = outputPrefix+".recovered_gtm";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	for(j = 0; j < num_muta; j++) {
		k = best_indices[j];
		for(i = 0; i < num_cell; i++) {
			if(i < num_cell-1) {
				ofs << states[k*num_cell+i] << " ";
			}
			else {
				ofs << states[k*num_cell+i] << endl;
			}
		}
	}
	ofs.close();
	
	fn = outputPrefix+".mutation_assignment";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	for(k = 0; k < num_cluster; k++) {
		for(j = 0; j < num_muta; j++) {
			if(best_indices[j] != k) {
				continue;
			}
			ofs << j << " ";
		}
		ofs << endl;
	}
	ofs.close();
	
	//char buf[500];
	//sprintf(buf, "%s_clusters_%d.scite_input", outputPrefix.c_str(), num_cluster);
	//fn = buf;
	fn = outputPrefix+".scite_input";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	
	for(i = 0; i < num_cluster; i++) {
		for(j = 0; j < num_cell; j++) {
			if(j < num_cell-1) {
				ofs << states[i*num_cell+j] << " ";
			}
			else {
				ofs << states[i*num_cell+j] << endl;
			}
		}
	}
	ofs.close();
	
}

void HetCaller::saveTree() {
	int i, j, k;
	ofstream ofs;
	string fn, outputPrefix = config.getStringPara("output");
	
	int num_cell = obsData.getROWS();
	int num_muta = obsData.getCOLS();
	int num_cluster = candi_s[best_s_indx].num_cluster;
	
	map<int, string> clone_labels;
	map<int, int> counts;
	int cells_per_line = max(5, (int) sqrt(num_cell/(num_cluster+1)));
	for(i = 0; i < num_cell; i++) {
		k = cell_assignments[i];
		if(clone_labels[k].empty()) {
			clone_labels[k] = cLabels[i];
			counts[k] = 1;
		}
		else {
			counts[k]++;
			if(counts[k] == cells_per_line) {
				counts[k] = 0;
				clone_labels[k] = clone_labels[k]+" "+cLabels[i]+"\\n";
			}
			else {
				if(counts[k] == 1) {
					clone_labels[k] = clone_labels[k]+cLabels[i];
				}
				else {
					clone_labels[k] = clone_labels[k]+" "+cLabels[i];
				}
			}
		}
	}
	
	int mutas_per_line = max(5, (int) sqrt(num_muta/num_cluster));
	map<int, string> node_labels;
	counts.clear();
	for(i = 0; i < num_muta; i++) {
		k = best_indices[i];
		if(node_labels[k].empty()) {
			node_labels[k] = mLabels[i];
			counts[k] = 1;
		}
		else {
			counts[k]++;
			if(counts[k] == mutas_per_line) {
				counts[k] = 0;
				node_labels[k] = node_labels[k]+" "+mLabels[i]+"\\n";
			}
			else {
				if(counts[k] == 1) {
					node_labels[k] = node_labels[k]+mLabels[i];
				}
				else {
					node_labels[k] = node_labels[k]+" "+mLabels[i];
				}
			}
		}
	}
	
	fn = outputPrefix+".dot";
	ofs.open(fn.c_str());
	if(!ofs.is_open()) {
		cerr << "Error: cannot open file " << fn << endl;
		exit(-1);
	}
	
	ofs << "digraph T {" << endl;
	ofs << "node [color=deeppink4, fontcolor=black, penwidth=2];" << endl;
	for(i = 0; i < num_cluster; i++) {
		ofs << i << " [label=\"" << node_labels[i] << "\"];" << endl;
	}
	ofs << i << " [label=\"\"];" << endl;
	
	ofs << "node [color=lightgrey, fontcolor=black, penwidth=2.5];" << endl;
	for(i = 0; i < num_cluster+1; i++) {
		if(!clone_labels[i].empty()) {
			ofs << i+num_cluster+1 << " [label=\"" << clone_labels[i] << "\"];" << endl;
		}
	}
	
	ofs << "edge [penwidth=1.5];" << endl;
	for(i = 0; i < num_cluster; i++) {
		ofs << clonal_tree[i] << "->" << i << ";" << endl;
	}
	
	for(i = 0; i < num_cluster+1; i++) {
		if(!clone_labels[i].empty()) {
			ofs << i << "->" << i+num_cluster+1 << ";" << endl;
		}
	}
	
	ofs << "}" << endl;
	
	
	ofs.close();
}


