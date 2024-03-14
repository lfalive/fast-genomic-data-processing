//
// Created by 梦想家xixi on 2021/10/13.
//

#include "ActivityProfileState.h"
#include <iostream>

double ActivityProfileState::isActiveProb() const { return activeProb; }

void ActivityProfileState::setIsActiveProb(const double aP) { this->activeProb = aP; }

ActivityProfileState::ActivityProfileState(SimpleInterval const &loc, const double activeProb)
		: loc(loc), activeProb(activeProb) {
	//Mutect2Utils::validateArg(this->loc.size() == 1, "Location for an ActivityProfileState must have to size 1 bp.");
	//Mutect2Utils::validateArg(resultValue >= 0, "Result value isn't null and its < 0, which is illegal");
}

ActivityProfileState::ActivityProfileState(const ActivityProfileState &activityProfileState)
		: loc(activityProfileState.loc), activeProb(activityProfileState.activeProb) {
	//Mutect2Utils::validateArg(loc.size() == 1, "Location for an ActivityProfileState must have to size 1 bp.");
	//Mutect2Utils::validateArg(resultValue >= 0, "Result value isn't null and its < 0, which is illegal");
}

ActivityProfileState::ActivityProfileState(const std::string &refName, hts_pos_t pos, double activeProb)
		: loc(refName, pos, pos), activeProb(activeProb) {
}

ActivityProfileState::ActivityProfileState(int refName, hts_pos_t pos, double activeProb)
		: loc(refName, pos, pos), activeProb(activeProb) {
}

int ActivityProfileState::getOffset(Locatable *regionStartLoc) {
	//Mutect2Utils::validateArg(regionStartLoc != nullptr, "Null object is not allowed here.");
	return loc.getStart() - regionStartLoc->getStart();
}

SimpleInterval &ActivityProfileState::getLoc() {
	return loc;
}

std::ostream &operator<<(std::ostream &os, ActivityProfileState &activityProfileState) {
	std::cout << "loc:  " << activityProfileState.getLoc() << std::endl << "activeProb:"
	          << activityProfileState.isActiveProb() << std::endl;
	return os;
}
