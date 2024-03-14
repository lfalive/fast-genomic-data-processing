//
// Created by cluster on 22-9-27.
//

#include "buildTreeUtils.h"
#include <iostream>
#include <deque>


trieNode *
buildTreeUtils::buildTreeWithHaplotype(const std::vector<std::shared_ptr<Haplotype>> &haplotypes, bool isFloat) {
	if (haplotypes.empty()) {
		return nullptr;
	}
	char *reference = nullptr;
	int referenceLength = 0;
	int record;
	for (int j = 0; j < haplotypes.size(); j++) {
		if (haplotypes[j]->getIsReference()) {
			reference = reinterpret_cast<char *>(haplotypes[j]->getBases().get());
			referenceLength = haplotypes[j]->getBasesLength();
			record = j;
			break;
		}
	}

	if (reference == nullptr) {
		throw std::invalid_argument("there is no refHaplotypes");
	}
	int avxLen = isFloat ? sizeof(float) * 8 : sizeof(double) * 8;
	int i = 0;
	auto *root = new trieNode();
	trieNode *father = root;
	auto *child = new trieNode({record});
	father->addChild(child);
	father = child;
	while (referenceLength - 1 > i * avxLen) {
		i++;
		child = new trieNode({record});
		father->addChild(child);
		father = child;
	}
	for (int j = 0; j < haplotypes.size(); j++) {
		if (haplotypes[j]->getIsReference()) {
			continue;
		}
		int baseLen = haplotypes[j]->getBasesLength();
		char *bases = reinterpret_cast<char *>(haplotypes[j]->getBases().get());
		i = 0;
		father = root;
		bool flag = false;
		for (const auto &node: father->getChild()) {
			char *nodebases = reinterpret_cast<char *>(haplotypes[node->getIndex()[0]]->getBases().get());
			if (nodebases[0] == bases[0]) {
				father = node;
				node->addIndex(j);
				flag = true;
				break;
			}
		}
		if (!flag) {
			child = new trieNode({j});
			father->addChild(child);
			father = child;
		}
		while (1 + i * avxLen < baseLen - avxLen) {
			i++;
			bool flag = false;
			trieNode *tmp;
			for (const auto &node: father->getChild()) {
				if (i * avxLen < baseLen) {
					char *nodebases = reinterpret_cast<char *>(haplotypes[node->getIndex()[0]]->getBases().get());
					int nodebaseLen = haplotypes[node->getIndex()[0]]->getBasesLength();
					if (i * avxLen < nodebaseLen &&
					    isEqual(nodebases + (i - 1) * avxLen + 1, bases + (i - 1) * avxLen + 1, avxLen + 1)) {
						flag = true;
						tmp = node;
						break;
					}
				} else {
					char *nodebases = reinterpret_cast<char *>(haplotypes[node->getIndex()[0]]->getBases().get());
					int nodebaseLen = haplotypes[node->getIndex()[0]]->getBasesLength();
					if (baseLen == nodebaseLen &&
					    isEqual(nodebases + (i - 1) * avxLen + 1, bases + (i - 1) * avxLen + 1,
					            baseLen - avxLen * (i - 1) - 1)) {
						throw std::invalid_argument("there are two same haplotypes");
					}
				}
			}
			if (!flag) {
				auto *child = new trieNode({j});
				father->addChild(child);
				father = child;
			} else {
				tmp->addIndex(j);
				father = tmp;
			}
		}
		for (const auto &node: father->getChild()) {
			if (1 + i * avxLen < baseLen) {
				char *nodebases = reinterpret_cast<char *>(haplotypes[node->getIndex()[0]]->getBases().get());
				int nodebaseLen = haplotypes[node->getIndex()[0]]->getBasesLength();
				if (i * avxLen < nodebaseLen && isEqual(nodebases + i * avxLen + 1, bases + i * avxLen + 1, avxLen)) {
					break;
				}
			} else {
				char *nodebases = reinterpret_cast<char *>(haplotypes[node->getIndex()[0]]->getBases().get());
				int nodebaseLen = haplotypes[node->getIndex()[0]]->getBasesLength();
				if (baseLen == nodebaseLen &&
				    isEqual(nodebases + i * avxLen + 1, bases + i * avxLen + 1, baseLen - avxLen * i - 1)) {
					throw std::invalid_argument("there are two same haplotypes");
				}
			}
		}
		auto *child = new trieNode({j});
		father->addChild(child);
	}
	return root;
}

bool buildTreeUtils::isEqual(char *c1, char *c2, int len) {
	return memcmp(c1, c2, len) == 0;
}

void buildTreeUtils::deleteTree(trieNode *root) {
	if (root == nullptr)
		return;
	for (auto &node: root->getChild()) {
		deleteTree(node);
	}
	delete root;
}

void buildTreeUtils::printLayerTree(trieNode *root) {
	std::deque<trieNode *> records;
	records.emplace_back(root);
	while (!records.empty()) {
		std::deque<trieNode *> tmp;
		while (!records.empty()) {
			trieNode *node = records.front();
			records.pop_front();
			for (auto newnode: node->getChild()) {
				tmp.push_back(newnode);
				for (int i: newnode->getIndex()) {
					std::cout << i << ", ";
				}
				std::cout << "    ";
			}
		}
		std::cout << std::endl;
		records = tmp;
	}
}

trieNode *
buildTreeUtils::buildTreeWithHaplotype_same_height(const std::vector<std::shared_ptr<Haplotype>> &haplotypes,
                                                   bool isFloat) {
	if (haplotypes.empty()) {
		return nullptr;
	}
	char *reference = nullptr;
	int referenceLength = 0;
	int record;
	for (int j = 0; j < haplotypes.size(); j++) {
		if (haplotypes[j]->getIsReference()) {
			reference = reinterpret_cast<char *>(haplotypes[j]->getBases().get());
			referenceLength = haplotypes[j]->getBasesLength();
			record = j;
			break;
		}
	}

	if (reference == nullptr) {
		throw std::invalid_argument("there is no refHaplotypes");
	}
	int avxLen = isFloat ? sizeof(float) * 8 : sizeof(double) * 8;
	int i = 0;
	auto *root = new trieNode();
	trieNode *father = root;
	auto *child = new trieNode({record});
	father->addChild(child);
	father = child;
	while (referenceLength - 1 > i * avxLen) {
		i++;
		child = new trieNode({record});
		father->addChild(child);
		father = child;
	}
	for (int j = 0; j < haplotypes.size(); j++) {
		if (haplotypes[j]->getIsReference()) {
			continue;
		}
		int baseLen = haplotypes[j]->getBasesLength();
		char *bases = reinterpret_cast<char *>(haplotypes[j]->getBases().get());
		i = 0;
		father = root;
		bool flag = false;
		for (const auto &node: father->getChild()) {
			char *nodebases = reinterpret_cast<char *>(haplotypes[node->getIndex()[0]]->getBases().get());
			int nodebaseLen = haplotypes[node->getIndex()[0]]->getBasesLength();
			if (nodebases[0] == bases[0] && nodebaseLen == baseLen) {
				father = node;
				node->addIndex(j);
				flag = true;
				break;
			}
		}
		if (!flag) {
			child = new trieNode({j});
			father->addChild(child);
			father = child;
		}
		while (1 + i * avxLen < baseLen - avxLen) {
			i++;
			bool flag = false;
			trieNode *tmp;
			for (const auto &node: father->getChild()) {
				if (i * avxLen < baseLen) {
					char *nodebases = reinterpret_cast<char *>(haplotypes[node->getIndex()[0]]->getBases().get());
					int nodebaseLen = haplotypes[node->getIndex()[0]]->getBasesLength();
					if (i * avxLen < nodebaseLen &&
					    isEqual(nodebases + (i - 1) * avxLen + 1, bases + (i - 1) * avxLen + 1, avxLen + 1)) {
						flag = true;
						tmp = node;
						break;
					}
				} else {
					char *nodebases = reinterpret_cast<char *>(haplotypes[node->getIndex()[0]]->getBases().get());
					int nodebaseLen = haplotypes[node->getIndex()[0]]->getBasesLength();
					if (baseLen == nodebaseLen &&
					    isEqual(nodebases + (i - 1) * avxLen + 1, bases + (i - 1) * avxLen + 1,
					            baseLen - avxLen * (i - 1) - 1)) {
						throw std::invalid_argument("there are two same haplotypes");
					}
				}
			}
			if (!flag) {
				auto *child = new trieNode({j});
				father->addChild(child);
				father = child;
			} else {
				tmp->addIndex(j);
				father = tmp;
			}
		}
		for (const auto &node: father->getChild()) {
            char *nodebases = reinterpret_cast<char *>(haplotypes[node->getIndex()[0]]->getBases().get());
            int nodebaseLen = haplotypes[node->getIndex()[0]]->getBasesLength();
            if (baseLen == nodebaseLen &&
            isEqual(nodebases + i * avxLen + 1, bases + i * avxLen + 1, baseLen - avxLen * i - 1)) {
                throw std::invalid_argument("there are two same haplotypes");
            }
		}
		auto *child = new trieNode({j});
		father->addChild(child);
	}
	return root;
}

unsigned long long buildTreeUtils::numberOfNodes(trieNode *root) {
	if (root == nullptr)
		return 0;
	unsigned long long ret = 1;
	for (const auto &node: root->getChild()) {
		ret += numberOfNodes(node);
	}
	return ret;
}
