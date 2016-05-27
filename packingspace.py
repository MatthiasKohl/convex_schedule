# TODO verify every function here with tests ! some are not good

# return whether or not the given coordinate lies inside the interval given by
# intervalCoord and intervalSize, according to the torus bound coordBound
# if isStrict, then the given coordinate must not lie on the boundaries of the interval
# all coordinates are assumed to be modulo coordBound
def isCoordIn(coord, intervalCoord, intervalSize, coordBound, isStrict=True):
    intervalCoord = intervalCoord % coordBound
    return any((isStrict and c > intervalCoord and c < intervalCoord + intervalSize)
               or (not isStrict and c >= intervalCoord and c <= intervalCoord + intervalSize)
               for c in [coord, coord + coordBound])

# coordinates represent the lowest possible coordinate in all dimensions of this space
# note that lowest possible does not necessarily mean the smallest number as we are in
# a torus topology
# boundaries represent the size of the space in all dimensions
# coordBoundaries represent the size of the topology in all dimensions
class ConvexSpace:
    def __init__(self, coordinates, boundaries, coordBoundaries):
        self.coordinates = coordinates
        self.boundaries = boundaries
        self.coordBoundaries = coordBoundaries

    def __hash__(self):
        return hash((tuple(self.coordBoundaries), tuple(self.coordinates), tuple(self.boundaries)))

    def __eq__(self, other):
        if isinstance(other, ConvexSpace):
            r = self.coordBoundaries == other.coordBoundaries
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

    # for all following functions, the coordinate bounds of the two spaces are assumed to be equa
    l
    # self contains other if all other coordinates lie inside the intervals defined by self's coords
    def contains(self, other):
        return all(isCoordIn(otherC, thisC, thisB, coordB) and
                   isCoordIn((otherC + otherB) % coordB, thisC, thisB, coordB)
                   for thisC, thisB, otherC, otherB, coordB in
                   zip(self.coordinates, self.boundaries, other.coordinates, other.boundaries,
                       self.coordBoundaries))

    # self intersects other if in all dimensions, one of the other coordinates lies in or on the
    # border of self's intervals
    def isIntersecting(self, other):
        # check if for all dimensions, one of the other coordinates lies inside our coordinates
        return all(isCoordIn(otherC, thisC, thisB, coordB, False) or
                   isCoordIn((otherC + otherB) % coordB, thisC, thisB, coordB, False)
                   for thisC, thisB, otherC, otherB, coordB in
                   zip(self.coordinates, self.boundaries, other.coordinates, other.boundaries,
                       self.coordBoundaries))

    # return intersection between this space and another as
    # a list of intervals in all dimensions or an empty list (if no intersection)
    def intersectionIntervals(self, other):
        intervals = []
        for thisC, thisB, otherC, otherB, coordB in
        zip(self.coordinates, self.boundaries, other.coordinates, other.boundaries,
            self.coordBoundaries):
            # align other interval such that this interval is at 0 (this way, we ignore thisC)
            otherCAligned = otherC - thisC if otherC >= thisC else otherC + coordB - thisC
            if (otherCAligned > thisB and otherCAligned + otherB < coordB):
                intervals = []
                break
            # un-align the coordinate
            newCoord = thisC if otherCAligned > thisB else (otherCAligned + thisC) % coordB
            # if otherCAligned is inside this interval, new size just goes until end of this interval
            # if not, it's the minimum of where the two intervals end
            # TODO this might not be correct
            if (thisB >= coordB):
                newSize = otherB
            else:
                newSize = min((otherCAligned + otherB) % coordB, thisB)
                if otherCAligned > thisB else thisB - otherCAligned
            intervals.append((newCoord, newSize))
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
            # align interval to ours
            if (interval[0] != self.coordinates[d]):
                newBoundaries = list(self.boundaries)
                newBoundaries[d] = interval[0] - self.coordinates[d] if interval[0] > self.coordinates[d]
                else interval[0] + self.coordBoundaries[d] - self.coordinates[d]
                newSpaces.append(ConvexSpace(list(self.coordinates), newBoundaries))
            if ((interval[0] + interval[1]) % self.coordBoundaries[d] !=
                (self.coordinates[d] + self.boundaries[d]) % self.coordBoundaries[d]):
                newCoordinates = list(self.coordinates)
                newCoordinates[d] = (interval[0] + interval[1]) % self.coordBoundaries[d]
                newBoundaries = list(self.boundaries)
                newBoundaries[d] = self.coordinates[d] + self.boundaries[d] - (interval[0] + interval[1])
                if self.coordinates[d] + self.boundaries[d] > interval[0] + interval[1] else
                self.coordinates[d] + self.boundaries[d] + self.coordBoundaries[d] - (interval[0] + interval[1])
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
        maxCoord = max(self.coordinates[i] + self.boundaries[i],
                       other.coordinates[i] + other.boundaries[i])
        joinedBoundaries[i] = maxCoord - joinedCoordinates[i] if maxCoord > joinedCoordinates[i]
        else maxCoord + self.coordBoundaries[i] - joinedCoordinates[i]
        # if joined interval fills entire space, align to 0
        if (joinedBoundaries[i] >= self.coordBoundaries[i]):
            joinedCoordinates[i] = 0
            joinedBoundaries[i] = self.coordBoundaries[i]
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
        maxCoord = max(self.coordinates[i] + self.boundaries[i],
                       other.coordinates[i] + other.boundaries[i])
        newBoundaries[i] = maxCoord - newCoordinates[i] if maxCoord > newCoordinates[i]
        else maxCoord + self.coordBoundaries[i] - newCoordinates[i]
        # if joined interval fills entire space, align to 0
        if (newBoundaries[i] >= self.coordBoundaries[i]):
            newCoordinates[i] = 0
            newBoundaries[i] = self.coordBoundaries[i]
        return ConvexSpace(newCoordinates, newBoundaries)

    def __str__(self):
        return str(self.coordinates) + ' -> ' + str(self.boundaries)
