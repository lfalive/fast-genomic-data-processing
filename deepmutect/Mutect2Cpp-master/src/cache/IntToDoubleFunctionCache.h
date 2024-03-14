//
// Created by 梦想家xixi on 2021/10/30.
//

#ifndef MUTECT2CPP_MASTER_INTTODOUBLEFUNCTIONCACHE_H
#define MUTECT2CPP_MASTER_INTTODOUBLEFUNCTIONCACHE_H


class IntToDoubleFunctionCache {
private:
	double *cache;
	int length;

protected:
	virtual int maxSize() = 0;

	virtual double compute(int n) = 0;

public:
	IntToDoubleFunctionCache() : length(0), cache(nullptr) {}

	double get(int i);

	void expandCache(int newCapacity);

	~IntToDoubleFunctionCache();
};


#endif //MUTECT2CPP_MASTER_INTTODOUBLEFUNCTIONCACHE_H
