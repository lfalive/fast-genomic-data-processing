//
// Created by hlf on 7/8/22.
//

#include "PairHMMConcurrentControl.h"

std::queue<std::shared_ptr<LikelihoodsTask>> PairHMMConcurrentControl::pairHMMTaskQueue;
std::queue<std::shared_ptr<LikelihoodsTask_trie>> PairHMMConcurrentControl::pairHMMTaskQueue_trie;
bool PairHMMConcurrentControl::startPairHMMConcurrentMode;
std::mutex PairHMMConcurrentControl::pairHMMMutex;
//std::atomic<unsigned long long> PairHMMConcurrentControl::trie_nodes;
//std::atomic<unsigned long long> PairHMMConcurrentControl::original_nodes;
//std::atomic<unsigned long long> PairHMMConcurrentControl::all_cases;
//std::atomic<unsigned long long> PairHMMConcurrentControl::unique_cases;
//std::atomic<unsigned long long> PairHMMConcurrentControl::all_reads;
//std::atomic<unsigned long long> PairHMMConcurrentControl::unique_reads;
//std::atomic<unsigned long long> PairHMMConcurrentControl::compute_double_cases;

void PairHMMConcurrentControl::initial() {
	pairHMMTaskQueue = std::queue<std::shared_ptr<LikelihoodsTask>>();
	pairHMMTaskQueue_trie = std::queue<std::shared_ptr<LikelihoodsTask_trie>>();
	startPairHMMConcurrentMode = false;
//	trie_nodes =  0;
//	original_nodes = 0;
//	all_cases = 0;
//	unique_cases = 0;
//	all_reads = 0;
//	unique_reads = 0;
//	compute_double_cases = 0;
}

LikelihoodsTask::LikelihoodsTask(std::vector<testcase> &taskTestcases, std::vector<double> &taskLikelihoodArray,
                                 unsigned long index, unsigned long testcasesSize)
		: taskTestcases(taskTestcases), taskLikelihoodArray(taskLikelihoodArray), index(index),
		  count(index), testcasesSize(testcasesSize) {}

LikelihoodsTask_trie::LikelihoodsTask_trie(std::vector<trie_testcase> &taskTestcases,
                                           std::vector<std::vector<double>> &taskLikelihoodArray, unsigned long index,
                                           unsigned long testcasesSize) : taskTestcases(taskTestcases),
                                                                          taskLikelihoodArray(taskLikelihoodArray),
                                                                          index(index), count(index),
                                                                          testcasesSize(testcasesSize) {}

