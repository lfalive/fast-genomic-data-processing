//
// Created by cluster on 22-9-27.
//

#ifndef MUTECT2CPP_MASTER_TRIENODE_H
#define MUTECT2CPP_MASTER_TRIENODE_H

#include <vector>

class trieNode {
private:
	int size;
	std::vector<int> index;
	std::vector<trieNode *> childs;

public:
	trieNode();

	~trieNode();

	trieNode(const std::vector<int> &index);

	void addChild(trieNode *node);

	std::vector<trieNode *> &getChild();

	trieNode &operator=(trieNode const &node);

	void addIndex(int child);

	std::vector<int> &getIndex();
};


#endif //MUTECT2CPP_MASTER_TRIENODE_H
