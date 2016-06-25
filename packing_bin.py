from packing_space import ConvexSpace
import itertools

# choose the best join for a flat fit
# simply try all possible joins and choose the one which leaves the biggest sizes
# according to reverse order dimensions, this is the biggest space in the order of
# dimensions
# returns the list of joined spaces
def best_flat_join(space1, space2):
    joined_spaces = []
    biggest_sizes = max(tuple(reversed(s.boundaries)) for s in [space1, space2])
    for d in range(len(space1.boundaries)):
        spaces = space1.join(space2, d)
        if (not spaces):
            continue
        max_size = max(tuple(reversed(s.boundaries)) for s in spaces)
        if (max_size >= biggest_sizes):
            # if a the biggest joined space has the same size than another, take it
            # (since the other could be the initial which should always be overridden if possible)
            biggest_sizes = max_size
            joined_spaces = spaces
    return joined_spaces

class Bin:
    def __init__(self, boundaries):
        self.boundaries = list(boundaries)
        self.freelist = set([ConvexSpace([0 for b in boundaries], list(boundaries),
                                         self.boundaries)])
        self.spaces = set()
        self.space_IDs = dict()

    def can_fit(self, boundaries_to_fit):
        return any(space.can_fit(boundaries_to_fit) for space in self.freelist)

    # fit the given boundaries in this bin over the i-th dimension (boundaries_to_fit should
    # be flat up to that dimension), then leave the biggest possible spaces in the order
    # of dimensions
    # (free spaces are sorted over all other dimensions first, then i-th dimension)
    # this assumes that can_fit was called before
    def fit_flat(self, boundaries_to_fit, i, ID = None):
        FS = min(filter(lambda s: s.can_fit(boundaries_to_fit), self.freelist),
                    key=lambda s:
                    [x for d, x in enumerate(s.coordinates) if d != i] + [s.coordinates[i]])
        AS = ConvexSpace(list(FS.coordinates), list(boundaries_to_fit), self.boundaries)
        #print('Fitting ' + str(boundaries_to_fit) + ' into ' + str(FS))
        self.spaces.add(AS)
        self.space_IDs[AS] = ID
        self.freelist.remove(FS)
        # print('Free space: ' + str(FS) + ' minus assigned ' + str(AS) +
        #       ': ' + str([str(x) for x in FS.minus(AS)]))
        new_free_spaces = set(FS.minus(AS))

        # join each pair of spaces together (according to best flat join)
        # until there is no update anymore (no pair of spaces could be joined)
        # TODO: this is not an efficient way of doing it since when a given space cannot
        # be joined to all other spaces, it doesn't need to be reconsidered anymore
        # a fast implementation would store a table with the different possible configurations
        # and do the entire loop (all joins) in one step based on that
        is_updated = True
        while(is_updated):
            is_updated = False
            for fs1, fs2 in itertools.combinations(sorted(new_free_spaces,
                                                          key=lambda s: tuple(reversed(s.boundaries)),
                                                          reverse = True), 2):
                #print('Trying to join ' + str(fs1) + ' and ' + str(fs2))
                joined_spaces = best_flat_join(fs1, fs2)
                if (joined_spaces):
                    #print('Joined spaces: ' + str([str(x) for x in joined_spaces]))
                    new_free_spaces.difference_update({fs1, fs2})
                    new_free_spaces.update(joined_spaces)
                    is_updated = True
                    break

        #print('New free spaces: ' + str([str(x) for x in new_free_spaces]))
        self.freelist.update(new_free_spaces)

    # fit the given boundaries at the best location according to caving degree
    # i.e. always fit at a corner and try to fill all dimensions as much as possible
    def fit_best(self, boundaries_to_fit, ID = None):
        # get the space with the minimum gap over the most dimensions
        # for each space, sort the boundaries by gap with boundariesToFIt non-decreasing
        # then take the minimum over that list (i.e. the space that has min minimal gap,
        # min 2nd minimal gap...)
        chosen_fs = min(filter(lambda s: s.can_fit(boundaries_to_fit), self.freelist),
                       key = lambda s:
                       sorted([b - boundaries_to_fit[d] for d, b in enumerate(s.boundaries)]))
        # TODO coordinates could be adapted depending on other free spaces
        # in order to leave more or less convex free space open
        chosen_fs_coords = list(chosen_fs.coordinates)
        AS = ConvexSpace(chosen_fs_coords, list(boundaries_to_fit), self.boundaries)
        self.spaces.add(AS)
        self.space_IDs[AS] = ID

        #print('Fitted bin is as follows:\n' + str(self))
        # make sure that the assigned space is subtracted from all existing FSs
        removed_fs = set()
        added_fs = set()
        for space in self.freelist:
            if (space.is_intersecting(AS)):
                removed_fs.add(space)
                added_fs.update(space.minus(AS))
        self.freelist.difference_update(removed_fs)
        self.freelist.update(added_fs)

        #print('Fitted bin is as follows:\n' + str(self))
        # make sure that any adjacent FSs are joined until no adjacent FSs can be joined
        # at most this can happen as much as there are FSs in freelist
        while(True):
            #print('Fitted bin is as follows:\n' + str(self))
            joined_space = None
            for fs1, fs2 in itertools.combinations(self.freelist, 2):
                joined_space = fs1.join_adjacent(fs2)
                if (joined_space != None):
                    break
            if (joined_space == None):
                break
            self.freelist.difference_update({fs1, fs2})
            self.freelist.add(joined_space)

        #print('Fitted space ' + str(boundaries_to_fit) + ' into bin as follows:\n' + str(self))

    # simple test to check if a bin configuration is possible according to given boundaries
    # it can be specified that free spaces must not overlap (usually, they can overlap
    # between each other, not with other spaces)
    def test_possible(self, allow_overlapping_fs = True, allow_adjacent_fs = True):
        # check that no space goes outside of given boundaries
        all_spaces = list(self.freelist) + list(self.spaces)
        BS = ConvexSpace([0 for x in self.boundaries], self.boundaries, self.boundaries)
        if (any(not BS.contains(s) for s in all_spaces)):
            print('Space ' + str([s for s in all_spaces if not BS.contains(s)][0]) +
                  ' is out of bounds')
            return False
        # check that there is no intersection between any two spaces
        overlap_spaces = list(self.spaces) if allow_overlapping_fs else all_spaces
        for s1, s2 in itertools.combinations(overlap_spaces, 2):
            if (s1.is_intersecting(s2)):
                print('Space ' + str(s1) + ' is intersecting with ' + str(s2) + ' but must not')
                return False
        # if free spaces can overlap, still need to check that they don't overlap with assigned
        if (allow_overlapping_fs):
            for s1, s2 in itertools.product(self.freelist, self.spaces):
                if (s1.is_intersecting(s2)):
                    print('Space ' + str(s1) + ' is intersecting with ' + str(s2) + ' but must not')
                    return False
        if (not allow_adjacent_fs):
            # there must not be any two adjacent free spaces
            for s1, s2 in itertools.combinations(self.freelist, 2):
                if (s1.join_adjacent(s2) != None):
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
