//
// Created by cluster on 4/3/23.
//

#ifndef MUTECT2CPP_MASTER_GRAPHOBJECTPOOL_H
#define MUTECT2CPP_MASTER_GRAPHOBJECTPOOL_H


#include <vector>
#include <memory>
#include <cassert>
#include "graph/MultiDeBruijnVertex.h"
#include "graph/MultiSampleEdge.h"
#include "graph/SeqVertex.h"
#include "graph/BaseEdge.h"

#define likely(x)    __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)

template<class T>
class ObjectPool {
private:
	size_t m_capacity;
	size_t m_capacity_initial;
	size_t m_capacity_inc;
	size_t m_cur;
	std::vector<T *> m_objects;

	void expand(size_t size) {
		try {
			m_capacity += size;
			m_objects.reserve(m_capacity);
			T *new_objects = new T[size];
			for (int i = 0; i < size; ++i)
				m_objects.emplace_back(&new_objects[i]);
		} catch (std::bad_alloc &) {
			throw std::runtime_error("ObjectPool: failed to expand");
		}
	}

public:
	ObjectPool(int initialCapacity, int inc) : m_capacity_initial(initialCapacity), m_capacity_inc(inc), m_capacity(0),
	                                           m_cur(0), m_objects(std::vector<T *>()) {
		assert(m_capacity_inc > 0);
		assert(m_capacity_initial > 0);
		expand(m_capacity_initial);
	}

	~ObjectPool() {
		delete[] m_objects[0];
		for (int block_p = m_capacity_initial; block_p < m_capacity; block_p += m_capacity_inc) {
			delete[] m_objects[block_p];
		}
		m_objects.clear();
//		std::cout << "~ObjectPool\n";
	}

	template<typename... Args>
	std::shared_ptr<T> create(Args &&... args) {
		if (unlikely(m_cur == m_capacity))
			expand(m_capacity_inc);

		try {
			m_objects[m_cur]->~T();
			new(m_objects[m_cur]) T(std::forward<Args>(args)...);
			return std::shared_ptr<T>(m_objects[m_cur++], [](T *p) {});
		} catch (...) {
			throw std::runtime_error("ObjectPool: failed to construct object");
		}
	}

	void reset() {
		m_cur = 0;
	}

	[[nodiscard]] size_t getCapacity() const {
		return m_capacity;
	}
};

class GraphObjectPool {
private:
	static thread_local ObjectPool<MultiDeBruijnVertex> MultiVertexPool;
	static thread_local ObjectPool<MultiSampleEdge> MultiEdgePool;
	static thread_local ObjectPool<SeqVertex> SeqVertexPool;
	static thread_local ObjectPool<BaseEdge> SeqEdgePool;
public:
	template<typename... Args>
	static std::shared_ptr<MultiDeBruijnVertex> createMultiVertex(Args &&... args) {
		return MultiVertexPool.create(std::forward<Args>(args)...);
	}

	template<typename... Args>
	static std::shared_ptr<MultiSampleEdge> createMultiEdge(Args &&... args) {
		return MultiEdgePool.create(std::forward<Args>(args)...);
	}

	template<typename... Args>
	static std::shared_ptr<SeqVertex> createSeqVertex(Args &&... args) {
		return SeqVertexPool.create(std::forward<Args>(args)...);
	}

	template<typename... Args>
	static std::shared_ptr<BaseEdge> createSeqEdge(Args &&... args) {
		return SeqEdgePool.create(std::forward<Args>(args)...);
	}

	static void reset(int thread_id);
};

#endif //MUTECT2CPP_MASTER_GRAPHOBJECTPOOL_H
