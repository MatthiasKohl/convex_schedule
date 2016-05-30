# return whether or not the interval contains the other interval according to bounds
def intervalContains(thisC, thisB, otherC, otherB, coordB, isStrict = False):
    if (thisB >= coordB):
        return True
    otherCAligned = otherC - thisC if otherC >= thisC else otherC + coordB - thisC
    if (isStrict):
        return otherCAligned != 0 and otherCAligned + otherB < thisB
    else:
        return otherCAligned + otherB <= thisB

# return whether or not the interval intersects the other interval according to bounds
def intervalIntersects(thisC, thisB, otherC, otherB, coordB, isStrict = True):
    if (otherB <= 0):
        return False
    otherCAligned = otherC - thisC if otherC >= thisC else otherC + coordB - thisC
    if (isStrict):
        return (otherCAligned < thisB) or (otherCAligned + otherB > coordB)
    else:
        return (otherCAligned <= thisB) or (otherCAligned + otherB >= coordB)

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
        return all(thisB >= otherB for thisB, otherB in zip(self.boundaries, boundariesToFit))

    # for all following functions, the coordinate bounds of the two spaces are assumed to be equal

    # self contains other if all other coordinates lie inside the intervals defined by self's coords
    def contains(self, other):
        return all(intervalContains(thisC, thisB, otherC, otherB, coordB)
                   for thisC, thisB, otherC, otherB, coordB in
                   zip(self.coordinates, self.boundaries, other.coordinates, other.boundaries,
                       self.coordBoundaries))

    # self intersects other if in all dimensions, one of the other coordinates lies in or on the
    # border of self's intervals
    def isIntersecting(self, other, isStrict=True):
        return all(intervalIntersects(thisC, thisB, otherC, otherB, coordB, isStrict)
                   for thisC, thisB, otherC, otherB, coordB in
                   zip(self.coordinates, self.boundaries, other.coordinates, other.boundaries,
                       self.coordBoundaries))

    # return intersection between this space and another as
    # a list of intervals in all dimensions or an empty list (if no intersection)
    # these intervals are relative to the coordinates of this space (since an interval can
    # wrap around this space). for instance, in 1-D, in a 24-sized torus,
    # the intervals (3, 21) and (23, 6) will have an intersection of (20, 3)
    def intersectionIntervals(self, other):
        intervals = []
        for thisC, thisB, otherC, otherB, coordB in\
        zip(self.coordinates, self.boundaries, other.coordinates, other.boundaries,
            self.coordBoundaries):
            # align other interval such that this interval is at 0 (this way, we ignore thisC)
            otherCAligned = otherC - thisC if otherC >= thisC else otherC + coordB - thisC
            if (otherCAligned > thisB and otherCAligned + otherB < coordB):
                intervals = []
                break
            # relative coordinate
            newCoord = 0 if otherCAligned >= thisB else otherCAligned
            # add sizes of the two possible intersections (end of this interval
            # and start of this interval)
            newSize = 0
            if (otherCAligned <= thisB):
                newSize = newSize + min(otherB, thisB - otherCAligned)
                if (otherCAligned + otherB > coordB):
                    newSize = newSize + (otherCAligned + otherB) % coordB
            else:
                newSize = min((otherCAligned + otherB) % coordB, thisB)
            intervals.append((newCoord, newSize))
        return intervals

    # in principle, each dimension can be cut in 2 new ConvexSpaces
    # this is unless there is no space left at any or both sides of this ConvexSpace
    # if there is no intersection (or empty intersection), then there are no new spaces
    def minus(self, other):
        newSpaces = []
        intervals = self.intersectionIntervals(other)
        if (any(s == 0 for c, s in intervals)):
            return newSpaces
        for d, interval in enumerate(intervals):
            coordB = self.coordBoundaries[d]
            # special case if this space is full torus interval
            if (self.boundaries[d] >= coordB and interval[1] < coordB):
                newCoordinates = list(self.coordinates)
                newCoordinates[d] = (self.coordinates[d] + interval[0] + interval[1]) % coordB
                newBoundaries = list(self.boundaries)
                newBoundaries[d] = coordB - interval[1]
                newSpaces.append(ConvexSpace(newCoordinates, newBoundaries, self.coordBoundaries))
                continue
            if (interval[0] != 0):
                newCoordinates = list(self.coordinates)
                newBoundaries = list(self.boundaries)
                # interval could wrap around, check for this case
                if (interval[0] + interval[1] > self.boundaries[d]):
                    newCoordinates[d] = ((interval[0] + interval[1]) % self.boundaries[d] +
                                         self.coordinates[d]) % coordB
                    newBoundaries[d] = newBoundaries[d] - interval[1]
                else:
                    newBoundaries[d] = interval[0]
                newSpaces.append(ConvexSpace(newCoordinates, newBoundaries, self.coordBoundaries))
            if (interval[0] + interval[1] < self.boundaries[d]):
                newCoordinates = list(self.coordinates)
                newCoordinates[d] = (self.coordinates[d] + interval[0] + interval[1]) % coordB
                newBoundaries = list(self.boundaries)
                newBoundaries[d] = self.boundaries[d] - (interval[0] + interval[1])
                newSpaces.append(ConvexSpace(newCoordinates, newBoundaries, self.coordBoundaries))
        return newSpaces

    # just get the joined space over i-th dimension according to intersection intervals
    # these intervals are assumed to not be empty and not 0 anywhere
    # (except for i-th dimension)
    def getJoinedSpace(self, other, i, intervals):
        # take full interval of both spaces in i-th dimension
        # and all other intervals as is from intersection to create joined space
        joinedCoordinates = [(x[0] + self.coordinates[d]) % self.coordBoundaries[d]
        for d, x in enumerate(intervals)]
        # the interval starts with the other coordinate if intersection starts at our coordinate
        if (intervals[i][0] == 0):
            joinedCoordinates[i] = other.coordinates[i]
        else:
            joinedCoordinates[i] = self.coordinates[i]

        joinedBoundaries = [x[1] for x in intervals]
        # the interval stops at other coordinate if intersection stops at our coordinate
        coordB = self.coordBoundaries[i]
        if (intervals[i][0] + intervals[i][1] == self.boundaries[i]):
            otherEnd = (other.coordinates[i] + other.boundaries[i]) % coordB
            alignment = 0 if otherEnd > joinedCoordinates[i] else coordB
            joinedBoundaries[i] = otherEnd + alignment - joinedCoordinates[i]
        else:
            selfEnd = (self.coordinates[i] + self.boundaries[i]) % coordB
            alignment = 0 if selfEnd > joinedCoordinates[i] else coordB
            joinedBoundaries[i] = selfEnd + alignment - joinedCoordinates[i]
        # if joined interval fills entire space, align to 0
        if (joinedBoundaries[i] >= coordB):
            joinedCoordinates[i] = 0
            joinedBoundaries[i] = coordB
        return ConvexSpace(joinedCoordinates, joinedBoundaries, self.coordBoundaries)

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

        joinedSpace = self.getJoinedSpace(other, i, intervals)
        # subtract the joinedSpace from the original spaces and return all resulting spaces
        return [joinedSpace] + self.minus(joinedSpace) + other.minus(joinedSpace)

    # return joined space containing self and other, if self is adjacent or overlapping in one
    # dimension and all other dimensions overlap entirely
    # return joined space if the update was possible, None if not
    def joinAdjacent(self, other):
        if (self.contains(other)):
            return self
        if (other.contains(self)):
            return other
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
        return self.getJoinedSpace(other, i, intervals)

    def __str__(self):
        return str(self.coordinates) + ' -> ' + str(self.boundaries)

def testSpaces():
    torusBounds = [8, 8]
    s1 = ConvexSpace([0, 0], [4, 4], torusBounds)
    s2 = ConvexSpace([2, 2], [4, 4], torusBounds)
    s3 = ConvexSpace([6, 6], [4, 4], torusBounds)
    s4 = ConvexSpace([2, 2], [6, 6], torusBounds)
    assert s4.contains(s2)
    assert not s2.contains(s4)
    assert not s4.contains(s3)
    assert not s3.contains(s4)
    assert not s2.isIntersecting(s3)
    assert all(x[1] == 0 for x in s2.intersectionIntervals(s3))
    assert s1.isIntersecting(s2)
    assert s1.intersectionIntervals(s2) == [(2, 2), (2, 2)]
    assert s1.intersectionIntervals(s3) == [(0, 2), (0, 2)]
    assert s1.getJoinedSpace(s2, 0, s1.intersectionIntervals(s2)) == ConvexSpace([0, 2], [6, 2], torusBounds)
    assert s2.getJoinedSpace(s1, 0, s2.intersectionIntervals(s1)) == ConvexSpace([0, 2], [6, 2], torusBounds)
    assert s1.intersectionIntervals(s4) == [(2, 2), (2, 2)]
    assert s2.intersectionIntervals(s3) == [(0, 0), (0, 0)]
    assert s3.intersectionIntervals(s1) == [(2, 2), (2, 2)]
    assert s1.getJoinedSpace(s3, 1, s1.intersectionIntervals(s3)) == ConvexSpace([0, 6], [2, 6], torusBounds)
    assert s3.getJoinedSpace(s1, 1, s3.intersectionIntervals(s1)) == ConvexSpace([0, 6], [2, 6], torusBounds)
