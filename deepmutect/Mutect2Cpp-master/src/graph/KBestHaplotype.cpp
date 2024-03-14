//
// Created by 梦想家xixi on 2021/11/24.
//

#include "KBestHaplotype.h"
#include <cmath>
#include <utility>

KBestHaplotype::KBestHaplotype(std::shared_ptr<SeqVertex> initialVertex,
                               std::shared_ptr<DirectedSpecifics<SeqVertex, BaseEdge>> graph)
		: Path<SeqVertex, BaseEdge>(std::move(initialVertex), std::move(graph)) {
	score = 0;
}

KBestHaplotype::KBestHaplotype(const std::shared_ptr<KBestHaplotype> &p, const std::shared_ptr<BaseEdge> &edge,
                               int totalOutgoingMultiplicity) : Path<SeqVertex, BaseEdge>(*p, edge) {
	score = p->getScore() + (double) (std::log10((long double) edge->getMultiplicity()) -
	                                  std::log10((long double) totalOutgoingMultiplicity));
	isReference &= edge->getIsRef();
}

std::shared_ptr<Haplotype> KBestHaplotype::getHaplotype() {
	int length = 0;
	std::shared_ptr<uint8_t[]> base = getBases(length);
	return std::make_shared<Haplotype>(base, length, getIsReference(), score);
}
