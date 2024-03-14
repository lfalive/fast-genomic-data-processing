//
// Created by cluster on 22-9-27.
//

#include <cstring>
#include "trieNode.h"
#include <iostream>

trieNode::trieNode() {
	size = 0;
}

trieNode &trieNode::operator=(trieNode const &node) {
	size = node.size;
	index = node.index;
	childs = node.childs;
	return *this;
}


std::vector<int> &trieNode::getIndex() {
	return index;
}

void trieNode::addIndex(int child) {
	index.emplace_back(child);
}

trieNode::~trieNode() {

}

trieNode::trieNode(const std::vector<int> &index) {
	this->index = index;
	size = index.size();
}

void trieNode::addChild(trieNode *node) {
	childs.emplace_back(node);
}

std::vector<trieNode *> &trieNode::getChild() {
	return childs;
}



