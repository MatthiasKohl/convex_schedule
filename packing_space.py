# return whether or not the first interval contains the second interval according to bounds
def interval_contains(c1, b1, c2, b2, cb, is_strict = False):
    if (b1 >= cb):
        return True
    c2_aligned = c2 - c1 if c2 >= c1 else c2 + cb - c1
    if (is_strict):
        return c2_aligned != 0 and c2_aligned + b2 < b1
    else:
        return c2_aligned + b2 <= b1

# return whether or not the given intervals intersect according to bounds
def intervalIntersects(c1, b1, c2, b2, cb, is_strict = True):
    if (b2 <= 0):
        return False
    c2_aligned = c2 - c1 if c2 >= c1 else c2 + cb - c1
    if (is_strict):
        return (c2_aligned < b1) or (c2_aligned + b2 > cb)
    else:
        return (c2_aligned <= b1) or (c2_aligned + b2 >= cb)

# coordinates represent the lowest possible coordinate in all dimensions of this space
# note that lowest possible does not necessarily mean the smallest number as we are in
# a torus topology
# boundaries represent the size of the space in all dimensions
# coord_boundaries represent the size of the topology in all dimensions
class ConvexSpace:
    def __init__(self, coordinates, boundaries, coord_boundaries):
        self.coordinates = coordinates
        self.boundaries = boundaries
        self.coord_boundaries = coord_boundaries

    def __hash__(self):
        return hash((tuple(self.coord_boundaries), tuple(self.coordinates), tuple(self.boundaries)))

    def __eq__(self, other):
        if isinstance(other, ConvexSpace):
            r = self.coord_boundaries == other.coord_boundaries
            r = r and self.coordinates == other.coordinates
            return r and self.boundaries == other.boundaries
        return NotImplemented

    def __ne__(self, other):
        result = self.__eq__(other)
        if result is NotImplemented:
            return result
        return not result

    def can_fit(self, boundaries):
        return not any(b1 < b2 for b1, b2 in zip(self.boundaries, boundaries))

    # for all following functions, the coordinate bounds of the two spaces are assumed to be equal

    # self contains other if all other coordinates lie inside the intervals defined by self's coords
    def contains(self, other):
        # same as all(interval_contains) but faster
        for c1, b1, c2, b2, cb in zip(self.coordinates, self.boundaries,
                                                        other.coordinates, other.boundaries,
                                                        self.coord_boundaries):
            if (b1 >= cb):
                continue
            if (c2 >= c1):
                if (c2 + b2 - c1 > b1):
                    return False
            elif (c2 + b2 + cb - c1 > b1):
                return False
        return True

    # self intersects other if in all dimensions, one of the other coordinates lies in or on the
    # border of self's intervals
    def is_intersecting(self, other, is_strict=True):
        # same as alll(intervalIntersects) but faster
        for c1, b1, c2, b2, cb in zip(self.coordinates, self.boundaries,
                                                        other.coordinates, other.boundaries,
                                                        self.coord_boundaries):
            if (c2 >= c1):
                if (c2 - c1 >= b1 and c2 - c1 + b2 <= cb):
                    return False
            elif (c2 + cb - c1 >= b1 and
                  c2 + cb - c1 + b2 <= cb):
                return False
        return True


    # return intersection between this space and another as
    # a list of intervals in all dimensions or an empty list (if no intersection)
    # these intervals are relative to the coordinates of this space (since an interval can
    # wrap around this space). for instance, in 1-D, in a 24-sized torus,
    # the intervals (3, 21) and (23, 6) will have an intersection of (20, 3)
    def intersection_intervals(self, other):
        intervals = []
        for c1, b1, c2, b2, cb in\
        zip(self.coordinates, self.boundaries, other.coordinates, other.boundaries,
            self.coord_boundaries):
            # align other interval such that this interval is at 0 (this way, we ignore c1)
            c2_aligned = c2 - c1 if c2 >= c1 else c2 + cb - c1
            if (c2_aligned > b1 and c2_aligned + b2 < cb):
                intervals = []
                break
            # relative coordinate (this is modulo the size of this interval)
            new_coord = 0 if c2_aligned >= b1 else c2_aligned
            # add sizes of the two possible intersections (end of this interval
            # and start of this interval)
            if (c2_aligned <= b1):
                new_boundary = min(b2, b1 - c2_aligned)
                if (c2_aligned + b2 > cb):
                    new_boundary = new_boundary + (c2_aligned + b2) % cb
            else:
                new_boundary = min((c2_aligned + b2) % cb, b1)
            intervals.append((new_coord, new_boundary))
        return intervals

    # in principle, each dimension can be cut in 2 new ConvexSpaces
    # this is unless there is no space left at any or both sides of this ConvexSpace
    # if there is no intersection (or empty intersection), then there are no new spaces
    def minus(self, other):
        new_spaces = []
        intervals = self.intersection_intervals(other)
        if (any(s == 0 for c, s in intervals)):
            return new_spaces
        for d, interval in enumerate(intervals):
            cb = self.coord_boundaries[d]
            # special case if this space is full torus interval
            if (self.boundaries[d] >= cb and interval[1] < cb):
                new_coordinates = list(self.coordinates)
                new_coordinates[d] = (self.coordinates[d] + interval[0] + interval[1]) % cb
                new_boundaries = list(self.boundaries)
                new_boundaries[d] = cb - interval[1]
                new_spaces.append(ConvexSpace(new_coordinates, new_boundaries, self.coord_boundaries))
                continue
            if (interval[0] != 0):
                new_coordinates = list(self.coordinates)
                new_boundaries = list(self.boundaries)
                # interval could wrap around, check for this case
                if (interval[0] + interval[1] > self.boundaries[d]):
                    new_coordinates[d] = ((interval[0] + interval[1]) % self.boundaries[d] +
                                         self.coordinates[d]) % cb
                    new_boundaries[d] = new_boundaries[d] - interval[1]
                else:
                    new_boundaries[d] = interval[0]
                new_spaces.append(ConvexSpace(new_coordinates, new_boundaries, self.coord_boundaries))
            if (interval[0] + interval[1] < self.boundaries[d]):
                new_coordinates = list(self.coordinates)
                new_coordinates[d] = (self.coordinates[d] + interval[0] + interval[1]) % cb
                new_boundaries = list(self.boundaries)
                new_boundaries[d] = self.boundaries[d] - (interval[0] + interval[1])
                new_spaces.append(ConvexSpace(new_coordinates, new_boundaries, self.coord_boundaries))
        return new_spaces

    # just get the joined space over i-th dimension according to intersection intervals
    # these intervals are assumed to not be empty and not 0 anywhere
    # (except for i-th dimension)
    def get_joined_space(self, other, i, intervals):
        cb = self.coord_boundaries[i]
        # take full interval of both spaces in i-th dimension
        # and all other intervals as is from intersection to create joined space
        joined_coords = [(x[0] + self.coordinates[d]) % self.coord_boundaries[d]
        for d, x in enumerate(intervals)]
        if (intervals[i][0] != 0 or (intervals[i][0] == 0 and intervals[i][1] == 0 and
            (self.coordinates[i] + self.boundaries[i]) % cb == other.coordinates[i])):
            # the interval starts with the other coordinate if intersection starts at our coordinate
            # special case when intervals touch at the right side of this interval
            joined_coords[i] = self.coordinates[i]
        else:
            joined_coords[i] = other.coordinates[i]

        joined_boundaries = [x[1] for x in intervals]
        # the interval stops at other coordinate if intersection stops at our coordinate
        # or in special case where intervals touch at right side of this interval
        if (intervals[i][0] + intervals[i][1] == self.boundaries[i] or
              (intervals[i][0] == 0 and intervals[i][1] == 0 and
               (self.coordinates[i] + self.boundaries[i]) % cb == other.coordinates[i])):
            other_end = (other.coordinates[i] + other.boundaries[i]) % cb
            alignment = 0 if other_end > joined_coords[i] else cb
            joined_boundaries[i] = other_end + alignment - joined_coords[i]
        else:
            self_end = (self.coordinates[i] + self.boundaries[i]) % cb
            alignment = 0 if self_end > joined_coords[i] else cb
            joined_boundaries[i] = self_end + alignment - joined_coords[i]
        # if joined interval fills entire space, align to 0
        if (joined_boundaries[i] >= cb):
            joined_coords[i] = 0
            joined_boundaries[i] = cb
        return ConvexSpace(joined_coords, joined_boundaries, self.coord_boundaries)

    # try to join two spaces over the i-th dimension
    def join(self, other, i):
        # if one space is contained in the other, return the greater one
        if (self.contains(other)):
            return [self]
        if (other.contains(self)):
            return [other]
        # check if all dimensions overlap (we can only join in that case)
        intervals = self.intersection_intervals(other)
        #print('Intersection intervals: ' + str(intervals))
        if (not intervals):
            return []
        # if any dimension of the intervals is 0 (except for the i-th),
        # return (can not join in that case)
        if (any([d != i and s == 0 for d, (c, s) in enumerate(intervals)])):
            return []

        joined_space = self.get_joined_space(other, i, intervals)
        # subtract the joined_space from the original spaces and return all resulting spaces
        return [joined_space] + self.minus(joined_space) + other.minus(joined_space)

    # return joined space containing self and other, if self is adjacent or overlapping in one
    # dimension and all other dimensions overlap entirely
    # return joined space if the update was possible, None if not
    def join_adjacent(self, other):
        intervals = self.intersection_intervals(other)
        if (not intervals):
            return None
        # get the dimension over which we join, all others have to overlap completely
        # at the same time, check if all intervals are contained in the same space
        i = -1
        contained = []
        for d in range(len(intervals)):
            if (intervals[d][1] != self.boundaries[d] or intervals[d][1] != other.boundaries[d]):
                contained.append(intervals[d][1] == self.boundaries[d])
                i = d
        # if only 1 interval was not overlapping, return the joined space
        # if all are contained in one of the spaces, return that one, else return None
        if (len(contained) == 1):
            return self.get_joined_space(other, i, intervals)
        if (all(contained)):
            return other
        if (any(contained)):
            return None
        return self

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
    assert not s2.is_intersecting(s3)
    assert all(x[1] == 0 for x in s2.intersection_intervals(s3))
    assert s1.is_intersecting(s2)
    assert s1.intersection_intervals(s2) == [(2, 2), (2, 2)]
    assert s1.intersection_intervals(s3) == [(0, 2), (0, 2)]
    assert s1.get_joined_space(s2, 0, s1.intersection_intervals(s2)) == ConvexSpace([0, 2], [6, 2], torusBounds)
    assert s2.get_joined_space(s1, 0, s2.intersection_intervals(s1)) == ConvexSpace([0, 2], [6, 2], torusBounds)
    assert s1.intersection_intervals(s4) == [(2, 2), (2, 2)]
    assert s2.intersection_intervals(s3) == [(0, 0), (0, 0)]
    assert s3.intersection_intervals(s1) == [(2, 2), (2, 2)]
    assert s1.get_joined_space(s3, 1, s1.intersection_intervals(s3)) == ConvexSpace([0, 6], [2, 6], torusBounds)
    assert s3.get_joined_space(s1, 1, s3.intersection_intervals(s1)) == ConvexSpace([0, 6], [2, 6], torusBounds)
