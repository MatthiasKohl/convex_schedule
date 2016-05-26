class ConvexSpace:
    def __init__(self, coordinates, boundaries, tag=0):
        self.coordinates = coordinates
        self.boundaries = boundaries
        self.tag = tag

    def __hash__(self):
        return hash((self.tag, tuple(self.coordinates), tuple(self.boundaries)))

    def __eq__(self, other):
        if isinstance(other, ConvexSpace):
            r = self.tag == other.tag
            r = r and self.coordinates == other.coordinates
            return r and self.boundaries == other.boundaries
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def canFit(self, boundariesToFit):
        for thisS, s in zip(self.boundaries, boundariesToFit):
            if (thisS < s):
                return False
        return True

    def contains(self, other):
        for thisCoord, thisD, otherCoord, otherD in zip(self.coordinates, self.boundaries,
                                                other.coordinates, other.boundaries):
            if (otherCoord < thisCoord or otherCoord > thisCoord + thisD
                or otherCoord + otherD > thisCoord + thisD):
                return False
        return True

    def isIntersecting(self, other):
        # check if for all dimensions, one of the other coordinates lies inside our coordinates
        for thisCoord, thisD, otherCoord, otherD in zip(self.coordinates, self.boundaries,
                                                        other.coordinates, other.boundaries):
            if (otherCoord + otherD <= thisCoord or otherCoord >= thisCoord + thisD):
                return False
        return True

    # return intersection between this space and another as
    # a list of intervals in all dimensions or an empty list (if no intersection)
    def intersectionIntervals(self, other):
        intervals = []
        for thisCoord, thisD, otherCoord, otherD in zip(self.coordinates, self.boundaries,
                                                        other.coordinates, other.boundaries):
            if (otherCoord + otherD < thisCoord or otherCoord > thisCoord + thisD):
                intervals = []
                break
            newCoord = max(thisCoord, otherCoord)
            newMaxCoord = min(thisCoord + thisD, otherCoord + otherD)
            intervals.append((newCoord, newMaxCoord - newCoord))
        return intervals

    # return the intersection of this space and another as a new space
    # if no intersection possible, the new space has dimensions 0
    def intersection(self, other):
        intervals = self.intersectionIntervals(other)
        if (not intervals):
            return ConvexSpace([0 for x in self.coordinates], [0 for x in self.boundaries])
        return ConvexSpace([x[0] for x in intervals], [x[1] for x in intervals])

    # in principle, each dimension can be cut in 2 new ConvexSpaces
    # this is unless there is no space left at any or both sides of this ConvexSpace
    # if there is no intersection (or empty intersection), then there are no new spaces
    def minus(self, other):
        newSpaces = []
        intervals = self.intersectionIntervals(other)
        if (any(s == 0 for c, s in intervals)):
            return newSpaces
        for d, interval in enumerate(intervals):
            if (interval[0] > self.coordinates[d]):
                newBoundaries = list(self.boundaries)
                newBoundaries[d] = interval[0] - self.coordinates[d]
                newSpaces.append(ConvexSpace(list(self.coordinates), newBoundaries))
            if (interval[0] + interval[1] < self.coordinates[d] + self.boundaries[d]):
                newCoordinates = list(self.coordinates)
                newCoordinates[d] = interval[0] + interval[1]
                newBoundaries = list(self.boundaries)
                newBoundaries[d] = self.coordinates[d] + self.boundaries[d] - (interval[0] + interval[1])
                newSpaces.append(ConvexSpace(newCoordinates, newBoundaries))
        return newSpaces

    # try to join two spaces over the i-th dimension
    def join(self, other, i):
        # if one space is contained in the other, return the greater one
        if (self.contains(other)):
            return [self]
        if (other.contains(self)):
            return [other]
        # check if all dimensions overlap (we can only join in that case)
        intervals = self.intersectionIntervals(other)
        #print('Intersection intervals: ' + str(intervals))
        if (not intervals):
            return []
        # if any dimension of the intervals is 0 (except for the i-th),
        # return (can not join in that case)
        if (any([d != i and s == 0 for d, (c, s) in enumerate(intervals)])):
            return []

        # take full interval of both spaces in i-th dimension
        # and all other intervals as is from intersection to create joined space
        joinedCoordinates = [x[0] for x in intervals]
        joinedCoordinates[i] = min(self.coordinates[i], other.coordinates[i])
        joinedBoundaries = [x[1] for x in intervals]
        joinedBoundaries[i] = max(self.coordinates[i] + self.boundaries[i],
                                  other.coordinates[i] + other.boundaries[i]) - joinedCoordinates[i]
        joinedSpace = ConvexSpace(joinedCoordinates, joinedBoundaries)

        # subtract the intersection from the original spaces and return all resulting spaces
        # which are not contained in the joined space already
        return [s for s in [joinedSpace] + self.minus(self.intersection(joinedSpace)) +
        other.minus(other.intersection(joinedSpace))
        if s == joinedSpace or not joinedSpace.contains(s)]

    # return joined space containing self and other, if self is adjacent or overlapping in one
    # dimension and all other dimensions overlap entirely
    # return joined space if the update was possible, None if not
    def joinAdjacent(self, other):
        if (self.contains(other)):
            return ConvexSpace(list(self.coordinates), list(self.boundaries))
        if (other.contains(self)):
            return ConvexSpace(list(other.coordinates), list(other.boundaries))
        intervals = self.intersectionIntervals(other)
        if (not intervals):
            return None
        # get the dimension over which we join, all others have to overlap completely
        i = -1
        for d in range(len(intervals)):
            if (intervals[d][1] != self.boundaries[d] or intervals[d][1] != other.boundaries[d]):
                if (i >= 0):
                    return None
                i = d
        newCoordinates = list(self.coordinates)
        newBoundaries = list(self.boundaries)
        newCoordinates[i] = min(self.coordinates[i], other.coordinates[i])
        newBoundaries[i] = max(self.coordinates[i] + self.boundaries[i],
                                      other.coordinates[i] + other.boundaries[i]) - newCoordinates[i]
        return ConvexSpace(newCoordinates, newBoundaries)

    def __str__(self):
        return str(self.coordinates) + ' -> ' + str(self.boundaries)
