//
// Created by cluster on 22-9-27.
//

#ifndef MUTECT2CPP_MASTER_BUILDTREEUTILS_H
#define MUTECT2CPP_MASTER_BUILDTREEUTILS_H

#include "Haplotype.h"
#include "trieNode.h"
#include <vector>

class buildTreeUtils {
public:
	static trieNode *buildTreeWithHaplotype(const std::vector<std::shared_ptr<Haplotype>> &haplotypes, bool isFloat);

	static trieNode *
	buildTreeWithHaplotype_same_height(const std::vector<std::shared_ptr<Haplotype>> &haplotypes, bool isFloat);

	static void printLayerTree(trieNode *root);

	static void deleteTree(trieNode *root);

	static unsigned long long numberOfNodes(trieNode *root);

private:
	static bool isEqual(char *c1, char *c2, int len);
};


#endif //MUTECT2CPP_MASTER_BUILDTREEUTILS_H
