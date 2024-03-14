//
// Created by cluster on 4/8/23.
//

#include "GraphObjectPool.h"

thread_local ObjectPool<MultiDeBruijnVertex> GraphObjectPool::MultiVertexPool(10000, 5000);
thread_local ObjectPool<MultiSampleEdge> GraphObjectPool::MultiEdgePool(10000, 5000);
thread_local ObjectPool<SeqVertex> GraphObjectPool::SeqVertexPool(10000, 5000);
thread_local ObjectPool<BaseEdge> GraphObjectPool::SeqEdgePool(10000, 5000);

void GraphObjectPool::reset(int thread_id) {
//	std::cout << "================== " << thread_id << " ==================\n";
//	std::cout << "MultiVertexPool\tcapacity: " << MultiVertexPool.getCapacity() << std::endl;
//	std::cout << "MultiEdgePool\tcapacity: " << MultiEdgePool.getCapacity() << std::endl;
//	std::cout << "SeqVertexPool\tcapacity: " << SeqVertexPool.getCapacity() << std::endl;
//	std::cout << "SeqEdgePool\tcapacity: " << SeqEdgePool.getCapacity() << std::endl;
	MultiVertexPool.reset();
	MultiEdgePool.reset();
	SeqVertexPool.reset();
	SeqEdgePool.reset();
//	std::cout << "GraphObjectPool clear.\n" << "================== " << thread_id << " ==================\n";
}
