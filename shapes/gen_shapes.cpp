#include <iostream>
#include <vector>
#include <deque>
#include <sstream>
#include <algorithm>

void prime_factor_table(int n, int *table) {
    if (n < 1) return;
    table[1] = 1;
    int i;
    for (i=2; i <= n; i++) {
        // trial division to fill the table
        int factor = 2;
        while(factor <= i) {
            if (i % factor == 0) {
                table[i] = factor;
                break;
            }
            factor++;
        }
    }
}

void prime_factors_rec(int* primeTable, int n, std::vector<int> *primes) {
    if (n <= 1) return;
    int factor = primeTable[n];
    (*primes).push_back(factor);
    prime_factors_rec(primeTable, n / factor, primes);
}

std::vector<int> prime_factors(int* primeTable, int n) {
    std::vector<int> primes;
    prime_factors_rec(primeTable, n, &primes);
    return primes;
}

bool inc_assignment(std::vector<int> *assignment, int nDim) {
    if ((*assignment).size() <= 0) return false;
    (*assignment)[0]++;
    size_t overflow = 0;
    while (overflow < (*assignment).size() && (*assignment)[overflow] >= nDim) {
        (*assignment)[overflow] = 0;
        overflow++;
        if (overflow < (*assignment).size()) (*assignment)[overflow]++;
    }
    return overflow < (*assignment).size();
}

template <class Deque, class T> void insert_into_deque(Deque& v, const T& t) {
  typename Deque::iterator i = std::lower_bound(v.begin(), v.end(), t);
  if (i == v.end() || t < *i)
    v.insert(i, t);
}

std::deque< std::vector<int> > shape_candidates(bool isGrid, int nDim, int *boundaries, int n, int *primeTable) {
    std::vector<int> primes = prime_factors(primeTable, n);
    std::deque< std::vector<int> > candidates;
    std::vector<int>::iterator primesMax = std::max_element(primes.begin(), primes.end());
    int *sizeMax = std::max_element(boundaries, boundaries+nDim);
    if (primesMax != primes.end() && sizeMax != boundaries + nDim && *primesMax > *sizeMax)
        return candidates;
    std::vector<int> assignment(primes.size(), 0);
    int convexDimSizeDivisor = isGrid ? 1 : 2;
    do {
        std::vector<int> distribution(nDim, 1);
        bool isOutOfRange = false;
        for (size_t i = 0; i < primes.size(); i++) {
            distribution[assignment[i]] *= primes[i];
            // a size equal to the size of the dimension is never out of range
            // for any other size, it must not be greater than the size of the dimension / divisor for convexity
            if (distribution[assignment[i]] != boundaries[assignment[i]] && distribution[assignment[i]] > boundaries[assignment[i]] / convexDimSizeDivisor) {
                isOutOfRange = true;
                break;
            }
        }
        if (!isOutOfRange)
            insert_into_deque(candidates, distribution);
    } while (inc_assignment(&assignment, nDim));
    return candidates;
}

/**
* Generate convex shapes for a grid or torus according to the following arguments:
* Arg 1: '0' if shapes are to be generated for a torus, '1' for a grid
* Arg 2-n: boundaries in all dimensions of the torus or grid as positive integers.
* There must be exactly n-1 dimensions with n > 1
*
* The generated shapes are written to stdout if everything is successful
* They are written in the CSV format, with n columns: the first column
* specifies the size of the shape, the n-1 following columns the boundaries
* in each dimension (same order of dimensions as in arguments 2-n)
*/
int main(int argc, char *argv[]) {
    if (argc < 3) return 2;
    int sizesOffset = 2;
    std::istringstream ssTop(argv[1]);
    bool isGrid;
    if (!(ssTop >> isGrid))
        return 2;
    int nDim = argc-sizesOffset;
    int boundaries[nDim];
    for (int i = 0; i < nDim; i++) {
        std::istringstream ss(argv[i+sizesOffset]);
        int x;
        if (!(ss >> x) || x < 1)
            return 2;
        boundaries[i] = x;
    }
    int total = 1;
    for (int i = 0; i < nDim; i++)
        total *= boundaries[i];
    int primeTable[total+1];
    prime_factor_table(total, primeTable);

    for (int i = 1; i <= total; i++) {
        std::deque< std::vector<int> > candidates =
            shape_candidates(isGrid, nDim, boundaries, i, primeTable);
        //std::cout << i << " -> " << candidates.size() << "\n";
        for (size_t j = 0; j < candidates.size(); j++) {
            std::cout << i << ",";
            for (size_t k = 0; k < candidates[j].size(); k++) {
                std::cout << candidates[j][k];
                if (k < candidates[j].size() - 1)
                    std::cout << ",";
            }
            std::cout << "\n";
        }
    }
    return 0;
}
