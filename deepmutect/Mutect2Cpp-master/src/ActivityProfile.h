//
// Class holding information about per-base activity scores for
// assembly region traversal
// Created by lhh on 10/25/21.
//

#ifndef MUTECT2CPP_MASTER_ACTIVITYPROFILE_H
#define MUTECT2CPP_MASTER_ACTIVITYPROFILE_H

#include <vector>
#include <optional>
#include <memory>
#include "htslib/sam.h"
#include "ActivityProfileState.h"
#include "AssemblyRegion.h"
#include "samtools/SAMFileHeader.h"

using namespace std;

class ActivityProfile {
private:
	/**
	 * Incorporate a single activity profile state into the current list of states
	 *
	 * If state's position occurs immediately after the last position in this profile, then
	 * the state is appended to the state list. If it's within the existing states list,
	 * the prob of stateToAdd is added to its corresponding state in the list. If the
	 * position would be the start of this profile, stateToAdd is simply ignored.
	 *
	 * @param stateToAdd
	 */
	void incorporateSingleState(const std::shared_ptr<ActivityProfileState> &stateToAdd);

	//int count = 0;

protected:
	deque<std::shared_ptr<ActivityProfileState>> stateList;
	int maxProbPropagationDistance;
	double activeProbThreshold;
	SimpleInterval regionStartLoc;
	SimpleInterval regionStopLoc;
	SAMFileHeader *header;

	int contigLength;

	vector<std::shared_ptr<AssemblyRegion>> *regions;

	optional<SimpleInterval> getLocForOffset(const SimpleInterval &relativeLoc, int offset);

	int findEndOfRegion(bool isActiveRegion, int minRegionSize, int maxRegionSize, bool forceConversion);

	/**
	 * Find the first index into the state list where the state is considered ! isActiveRegion
	 *
	 * Note that each state has a probability of being active, and this function thresholds that
	 * value on activeProbThreshold, coloring each state as active or inactive. Finds the
	 * largest contiguous stretch of states starting at the first state (index 0) with the same isActive
	 * state as isActiveRegion. If the entire state list has the same isActive value, then returns maxRegionSize
	 *
	 * @param isActiveRegion
	 * @param maxRegionSize
	 * @return
	 */
	int findFirstActivityBoundary(bool isActiveRegion, int maxRegionSize);

	/**
	 * Find the local minimum within 0 - endOfActiveRegion where we should divide region
	 *
	 * This algorithm finds the global minimum probability state within the region [minRegionSize, endOfActiveRegion)
	 * (exclusive of endOfActiveRegion), and returns the state index of that state
	 *
	 * @param endOfActiveRegion
	 * @param minRegionSize
	 * @return
	 */
	int findBestCutSite(int endOfActiveRegion, int minRegionSize);

	/**
	 * Is the probability at index in a local minimum?
	 * @param index
	 * @return
	 */
	bool isMinimum(int index);

public:
	ActivityProfile(int maxProbPropagationDistance, double activeProbThreshold, SAMFileHeader *header);

	// A virtual destructor is necessary for an abstract class. Otherwise delete operation will invoke undefined behavior
	virtual ~ActivityProfile();

	/**
	 * Is this profile empty?
	 * @return true if the profile is empty
	 * */
	bool isEmpty();

	/**
	 * Add the next ActivityProfileState to this profile
	 * @param state a well-formed ActivityProfileState result to incorporate into this profile
	 */
	void add(const std::shared_ptr<ActivityProfileState> &state);

	virtual vector<std::shared_ptr<ActivityProfileState>> *
	processState(const std::shared_ptr<ActivityProfileState> &justAddedState);

	int getEnd();

	/**
	 * Get the next completed assembly regions from this profile, and remove all states supporting them from this profile
	 *
	 * Takes the current profile and finds all of the active / inactive from the start of the profile that are
	 * ready. By ready we mean unable to have their probability modified any longer by future additions to the
	 * profile. The regions that are popped off the profile take their states with them, so the start of this
	 * profile will always be after the end of the last region returned here.
	 *
	 * This function may not return anything, if no regions are ready
	 *
	 * No returned regions will be larger than maxRegionSize.
	 * @param assemblyRegionExtension
	 * @param minRegionSize
	 * @param maxRegionSize
	 * @param forceConversion
	 * @return
	 */
	vector<std::shared_ptr<AssemblyRegion>> *
	popReadyAssemblyRegions(int assemblyRegionExtension, int minRegionSize, int maxRegionSize, bool forceConversion);

	std::shared_ptr<AssemblyRegion>
	popReadyAssemblyRegion(int assemblyRegionExtension, int minRegionSize, int maxRegionSize, bool forceConversion);

	virtual int getMaxProbPropagationDistance() = 0;

	void clear();
};


#endif //MUTECT2CPP_MASTER_ACTIVITYPROFILE_H
