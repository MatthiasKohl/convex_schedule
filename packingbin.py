from packingspace import ConvexSpace
import itertools

# choose the best join for a flat fit
# simply try all possible joins and choose the one which leaves the biggest sizes
# according to reverse order dimensions, this is the biggest space in the order of
# dimensions
# returns the list of joined spaces
def chooseBestJoinFlat(space1, space2):
    joinedSpaces = []
    biggestSizes = max(tuple(reversed(s.boundaries)) for s in [space1, space2])
    for d in range(len(space1.boundaries)):
        spaces = space1.join(space2, d)
        if (not spaces):
            continue
        maxSize = max(tuple(reversed(s.boundaries)) for s in spaces)
        if (maxSize >= biggestSizes):
            # if a the biggest joined space has the same size than another, take it
            # (since the other could be the initial which should always be overridden if possible)
            biggestSizes = maxSize
            joinedSpaces = spaces
    return joinedSpaces

class Bin:
    def __init__(self, boundaries):
        self.boundaries = boundaries
        self.freelist = set([ConvexSpace([0 for d in boundaries], list(boundaries))])
        self.spaces = set()

    def canFit(self, boundariesToFit):
        return any(space.canFit(boundariesToFit) for space in self.freelist)

    # fit the given boundaries in this bin over the i-th dimension (boundariesToFit should
    # be flat up to that dimension), then leave the biggest possible spaces in the order
    # of dimensions
    # (free spaces are sorted over all other dimensions first, then i-th dimension)
    # this assumes that canFit was called before
    def fitFlat(self, boundariesToFit, i):
        freeSpace = min(filter(lambda s: s.canFit(boundariesToFit), self.freelist),
                    key=lambda s: [x for d, x in enumerate(s.coordinates) if d != i] + [s.coordinates[i]])
        assignedSpace = ConvexSpace(list(freeSpace.coordinates), list(boundariesToFit))
        #print('Fitting ' + str(boundariesToFit) + ' into ' + str(freeSpace))
        self.spaces.add(assignedSpace)
        self.freelist.remove(freeSpace)
        newFreeSpaces = set(freeSpace.minus(assignedSpace))

        # join each pair of spaces together (according to best flat join)
        # until there is no update anymore (no pair of spaces could be joined)
        # TODO: this is not an efficient way of doing it since when a given space cannot
        # be joined to all other spaces, it doesn't need to be reconsidered anymore
        # a fast implementation would store a table with the different possible configurations
        # and do the entire loop (all joins) in one step based on that
        isUpdated = True
        while(isUpdated):
            isUpdated = False
            for fs1, fs2 in itertools.combinations(sorted(newFreeSpaces,
                                                          key=lambda s: tuple(reversed(s.boundaries)),
                                                          reverse = True), 2):
                #print('Trying to join ' + str(fs1) + ' and ' + str(fs2))
                joinedSpaces = chooseBestJoinFlat(fs1, fs2)
                if (joinedSpaces):
                    #print('Joined spaces: ' + str([str(x) for x in joinedSpaces]))
                    newFreeSpaces.difference_update({fs1, fs2})
                    newFreeSpaces.update(joinedSpaces)
                    isUpdated = True
                    break

        #print('New free spaces: ' + str([str(x) for x in newFreeSpaces]))
        self.freelist.update(newFreeSpaces)

    # fit the given boundaries at the best location according to caving degree
    # i.e. always fit at a corner and try to fill all dimensions as much as possible
    def fitBest(self, boundariesToFit):
        # get the space with the minimum gap over the most dimensions
        # for each space, sort the boundaries by gap with boundariesToFIt non-decreasing
        # then take the minimum over that list (i.e. the space that has min minimal gap,
        # min 2nd minimal gap...)
        chosenFS = min(filter(lambda s: s.canFit(boundariesToFit), self.freelist),
                       key = lambda s: sorted([b - boundariesToFit[d] for d, b in enumerate(s.boundaries)]))
        # TODO coordinates could be adapted depending on other free spaces
        # in order to leave more or less convex free space open
        chosenFSCoords = list(chosenFS.coordinates)
        assignedSpace = ConvexSpace(chosenFSCoords, list(boundariesToFit))
        self.spaces.add(assignedSpace)

        #print('Fitted bin is as follows:\n' + str(self))
        # make sure that the assigned space is subtracted from all existing FSs
        removedFSs = set()
        addedFSs = set()
        for space in self.freelist:
            if (space.isIntersecting(assignedSpace)):
                removedFSs.add(space)
                addedFSs.update(space.minus(assignedSpace))
        self.freelist.difference_update(removedFSs)
        self.freelist.update(addedFSs)

        #print('Fitted bin is as follows:\n' + str(self))
        # make sure that any adjacent FSs are joined until no adjacent FSs can be joined
        # at most this can happen as much as there are FSs in freelist
        while(True):
            #print('Fitted bin is as follows:\n' + str(self))
            removedFSs = set()
            addedFSs = set()
            for fs1, fs2 in itertools.combinations(self.freelist, 2):
                joinedSpace = fs1.joinAdjacent(fs2)
                if (joinedSpace != None):
                    removedFSs.update({fs1, fs2})
                    addedFSs.add(joinedSpace)
            if (not removedFSs and not addedFSs):
                break
            self.freelist.difference_update(removedFSs)
            self.freelist.update(addedFSs)

        #print('Fitted space ' + str(boundariesToFit) + ' into bin as follows:\n' + str(self))

    # simple test to check if a bin configuration is possible according to given boundaries
    # it can be specified that free spaces must not overlap (usually, they can overlap
    # between each other, not with other spaces)
    def testPossible(self, allowOverlappingFS = True, allowAdjFS = True):
        # check that no space goes outside of given boundaries
        allSpaces = list(self.freelist) + list(self.spaces)
        boundarySpace = ConvexSpace([0 for x in self.boundaries], self.boundaries)
        if (not all(boundarySpace.contains(s) for s in allSpaces)):
            print('Space ' + str([boundarySpace.contains(s) for s in allSpaces][0]) + ' is out of bounds')
            return False
        # check that there is no intersection between any two spaces
        overlapSpaces = list(self.spaces) if allowOverlappingFS else allSpaces
        for s1, s2 in itertools.combinations(overlapSpaces, 2):
            if (s1.isIntersecting(s2)):
                print('Space ' + str(s1) + ' is intersecting with ' + str(s2) + ' but must not')
                return False
        # if free spaces can overlap, still need to check that they don't overlap with assigned
        if (allowOverlappingFS):
            for s1, s2 in itertools.product(self.freelist, self.spaces):
                if (s1.isIntersecting(s2)):
                    print('Space ' + str(s1) + ' is intersecting with ' + str(s2) + ' but must not')
                    return False
        if (not allowAdjFS):
            # there must not be any two adjacent free spaces
            for s1, s2 in itertools.combinations(self.freelist, 2):
                if (s1.joinAdjacent(s2) != None):
                    print('Space ' + str(s1) + ' is adjacent with ' + str(s2) + ' but must not')
                    return False
        return True

    def __str__(self):
        s = ''
        for space in sorted(self.spaces, key=lambda s: s.coordinates):
            s = s + str(space) + '\n'
        s = s + 'Free list:\n'
        for space in sorted(self.freelist, key=lambda s: s.coordinates):
            s = s + str(space) + '\n'
        return s
