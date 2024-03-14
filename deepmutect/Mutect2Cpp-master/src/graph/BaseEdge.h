//
// Created by 梦想家xixi on 2021/10/20.
//

#ifndef MUTECT2CPP_MASTER_BASEEDGE_H
#define MUTECT2CPP_MASTER_BASEEDGE_H

#include <string>
#include "parallel_hashmap/phmap.h"
#include <memory>
#include <vector>

class BaseEdge {

private:
	int multiplicity;
	bool isRef;

public:
	BaseEdge(bool isRef, int multiplicity);

	BaseEdge() : isRef(false), multiplicity(0) {}

	virtual ~BaseEdge() = default;;

	/**
	 * Get the number of observations of paths connecting two vertices
	 * @return a positive integer >= 0
	 */
	int getMultiplicity() const { return multiplicity; }

	/**
	* Increase the multiplicity of this edge by incr
	* @param incr the change in this multiplicity, must be >= 0
	*/
	virtual void incMultiplicity(int incr);

	/**
	 * Set the multiplicity of this edge to value
	 * @param value an integer >= 0
	 */
	void setMultiplicity(int value);

	/**
	* Does this edge indicate a path through the reference graph?
	* @return true if so
	*/
	bool getIsRef() const { return isRef; }

	/**
	 * Indicate that this edge follows the reference sequence, or not
	 * @param ref true if this is a reference edge
	 */
	void setIsRef(bool ref);

	/**
	 * Add edge to this edge, updating isRef and multiplicity as appropriate
	 *
	 * isRef is simply the or of this and edge
	 * multiplicity is the sum
	 *
	 * @param edge the edge to add
	 * @return this
	 */
	BaseEdge add(BaseEdge &edge);

	static std::shared_ptr<BaseEdge> makeOREdge(const std::vector<std::shared_ptr<BaseEdge>> &edges, int multiplicity);

	static std::shared_ptr<BaseEdge>
	makeOREdge(const phmap::flat_hash_set<std::shared_ptr<BaseEdge>> &edges, int multiplicity);

};


#endif //MUTECT2CPP_MASTER_BASEEDGE_H
