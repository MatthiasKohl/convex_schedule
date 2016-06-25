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

void prime_factors_rec(int* prime_table, int n, std::vector<int> *primes) {
    if (n <= 1) return;
    int factor = prime_table[n];
    (*primes).push_back(factor);
    prime_factors_rec(prime_table, n / factor, primes);
}

std::vector<int> prime_factors(int* prime_table, int n) {
    std::vector<int> primes;
    prime_factors_rec(prime_table, n, &primes);
    return primes;
}

bool inc_assignment(std::vector<int> *assignment, int nDim) {
    if ((*assignment).size() <= 0) return false;
    (*assignment)[0]++;
    int overflow = 0;
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

std::deque< std::vector<int> > shape_candidates(bool isGrid, int nDim, int *dSizes, int n, int *prime_table) {
    std::vector<int> primes = prime_factors(prime_table, n);
    std::deque< std::vector<int> > candidates;
    std::vector<int>::iterator primesMax = std::max_element(primes.begin(), primes.end());
    int *sizeMax = std::max_element(dSizes, dSizes+nDim);
    if (primesMax != primes.end() && sizeMax != dSizes + nDim && *primesMax > *sizeMax)
        return candidates;
    std::vector<int> assignment(primes.size(), 0);
    int convexDimSizeDivisor = isGrid ? 1 : 2;
    do {
        std::vector<int> distribution(nDim, 1);
        bool isOutOfRange = false;
        for (int i = 0; i < primes.size(); i++) {
            distribution[assignment[i]] *= primes[i];
            // a size equal to the size of the dimension is never out of range
            // for any other size, it must not be greater than the size of the dimension / divisor for convexity
            if (distribution[assignment[i]] != dSizes[assignment[i]] && distribution[assignment[i]] > dSizes[assignment[i]] / convexDimSizeDivisor) {
                isOutOfRange = true;
                break;
            }
        }
        if (!isOutOfRange)
            insert_into_deque(candidates, distribution);
    } while (inc_assignment(&assignment, nDim));
    return candidates;
}

int main(int argc, char *argv[]) {
    if (argc < 3) return 2;
    int sizesOffset = 2;
    std::istringstream ssTop(argv[1]);
    bool isGrid;
    if (!(ssTop >> isGrid))
        return 2;
    int nDim = argc-sizesOffset;
    int dSizes[nDim];
    for (int i = 0; i < nDim; i++) {
        std::istringstream ss(argv[i+sizesOffset]);
        int x;
        if (!(ss >> x))
            return 2;
        dSizes[i] = x;
    }
    int total = 1;
    for (int i = 0; i < nDim; i++)
        total *= dSizes[i];
    int prime_table[total+1];
    prime_factor_table(total, prime_table);

    for (int i = 1; i <= total; i++) {
        std::deque< std::vector<int> > candidates =
            shape_candidates(isGrid, nDim, dSizes, i, prime_table);
        //std::cout << i << " -> " << candidates.size() << "\n";
        for (int j = 0; j < candidates.size(); j++) {
            std::cout << i << ",";
            for (int k = 0; k < candidates[j].size(); k++) {
                std::cout << candidates[j][k];
                if (k < candidates[j].size() - 1)
                    std::cout << ",";
            }
            std::cout << "\n";
        }
    }
    return 0;
}
