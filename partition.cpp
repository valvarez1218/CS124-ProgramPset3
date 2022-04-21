#include <iostream>
#include <iomanip>
#include <vector>
// #include <utility>
#include <fstream>
#include <time.h>
#include <cstdlib>
#include <assert.h>
#include <random>
#include <climits>
#include <cmath>

using namespace std;

int N = 100;
int g_NumIterations = 25000;

// structure for min heap
struct MaxHeap {

private:
    size_t heapSize;

public:
    MaxHeap();
    // vector representing heap
    vector<long> heap;

    // swap location of elements at indices i and j
    void swap(int, int);

    // given an index, sift element at index up
    void siftUp(int);

    // when we call getMin() we add last element to top and sift down
    void siftDown(int);
    
    // add element to end of heap and sift up
    void insert(long);

    // remove smallest element from heap
    long getMax();

    // return size of heap
    size_t size();

    // return whether or not the heap is empty
    bool empty();

    // clear the heap
    void clear();
};

// Constructor function
MaxHeap::MaxHeap() {
    // initialize heap size to 0
    heapSize = 0;
}

long KarmarkarKarp(MaxHeap &A) {
    while (true) {
        long val1 = A.getMax();
        
        // if the heap is empty after removing first value
        //      then this value is the residue
        if (A.empty()) {
            return val1;
        }

        // otherwise retrieve next value as second largest number
        long val2 = A.getMax();

        long difference = val1 - val2;

        // insert the difference into the heap
        A.insert(difference);
    }
}

// calculate residue when partition is represented as list S of 1 and -1
long plusMinusSum(vector<int> S, MaxHeap &A) {
    long residue = 0;
    assert(S.size() == A.size());

    for (int i=0; i < int(S.size()); i++) {
        residue += A.heap[i]*S[i];
    }
    return abs(residue);
}

vector<int> randMovePM(vector<int> &S, int numChanges) {
    vector<int> S_p = S;
    // srand(time(NULL));
    for (int k = 0; k < numChanges; k++) {
        int i = rand() % N;
        S_p[i] = -1 * S_p[i];
    }

    return S_p;
}

vector<int> randMovePrepartition(vector<int> &P, int numChanges) {
    vector<int> P_new = P;
    // srand(time(NULL));
    for (int k = 0; k < numChanges; k++) {
        int i = rand() % N;
        // generate j such that j != i
        int j = rand() % N;
        while (P_new[i] == j) {
            i = rand() % N;
        }
        // set p_i = j
        P_new[i] = j;
    }

    return P_new;
}

// given a prepartition, generate A' using algorithm as described in pset
MaxHeap generateA_p(MaxHeap &A, vector<int> &P) {
    vector<long> V_p(N, 0);
    assert(A.size() == P.size());

    for (int j = 0; j < N; j++) {
        int p_j = P[j];
        V_p[p_j] = V_p[p_j] + A.heap[j];
    }

    MaxHeap A_p;
    for (long entry : V_p) {
        A_p.insert(entry);
    }

    return A_p;
}

// generate N random numbers in range [1, 10^12]
void populateA(MaxHeap &A) {
    default_random_engine generator(random_device{}());
    long MAX = pow(10,12);
    uniform_int_distribution<long> distribution(1,MAX);
    for (int i=0; i < N; i++) {
        long r = distribution(generator);
        A.insert(r);
    }
}

// generate random sequence of length N of 1 and -1
vector<int> genPlusMinus() {
    vector<int> S;
    // srand(time(NULL));
    for (int i=0; i < N; i++) {
        int r = rand() % 2 == 0 ? -1 : 1;
        S.push_back(r);
    }
    return S;
}

// generate random prepartition, list of length N of numbers in
//      range [0, N-1]
vector<int> genRandPartition() {
    vector<int> P;
    // srand(time(NULL));
    for (int i = 0; i < N; i++) {
        int r = rand() % N;
        P.push_back(r);
    }
    return P;
}

long repeatedRandom(MaxHeap &A, bool plusMinus) {
    long bestResidue = 0;
    if (plusMinus) {
        vector<int> S = genPlusMinus();
        bestResidue = plusMinusSum(S, A);

        for (int i = 0; i < g_NumIterations; i++) {
            vector<int> S_p = genPlusMinus();
            long residue = plusMinusSum(S_p, A);
            if (residue < bestResidue) {
                S = S_p;
                bestResidue = residue;
            }
        }
    } else {
        vector<int> P = genRandPartition();
        MaxHeap A_p = generateA_p(A, P);
        bestResidue = KarmarkarKarp(A_p);

        for (int i = 0; i < g_NumIterations; i++) {
            vector<int> P_r = genRandPartition();
            MaxHeap A_r = generateA_p(A, P_r);
            long residue = KarmarkarKarp(A_r);
            if (residue < bestResidue) {
                P = P_r;
                A_p = A_r;
                bestResidue = residue;
            }
        }
    }
    
    return bestResidue;
}


long hillClimbing(MaxHeap& A, bool plusMinus) {
    // list of 1 and -1 representation
    if (plusMinus) {
        vector<int> S = genPlusMinus();
        assert(S.size() == A.size());

        vector<int> bestS = S;
        long bestResidue = plusMinusSum(S, A);

        for (int i = 0; i < g_NumIterations; i++) {
            vector<int> S_p = randMovePM(S, 5);
            long currResidue = plusMinusSum(S_p, A);
            if (currResidue < bestResidue) {
                // only change S if we found a better solution
                S = S_p;
                bestResidue = currResidue;
                bestS = S_p;
            }
        }

        return bestResidue;
    }
    // prepartition representation
    else {
        vector<int> P = genRandPartition();
        assert(P.size() == A.size());
        
        vector<int> bestP = P;
        MaxHeap A_p = generateA_p(A, P);
        long bestResidue = KarmarkarKarp(A_p);

        for (int i = 0; i < g_NumIterations; i++) {
            vector<int> P_new = randMovePrepartition(P, 5);
            A_p = generateA_p(A, P_new);
            long currResidue = KarmarkarKarp(A_p);

            if (currResidue < bestResidue) {
                // only change P if we found a better solution
                P = P_new;
                bestResidue = currResidue;
                bestP = P_new;
            }
        }

        return bestResidue;
    }
}

double prob_eT(int i, long resDiff) {
    double power = resDiff / pow(10, 10)*pow(0.8, floor(i/300.));

    return exp(power);
}

long simulatedAnnealing(MaxHeap& A, bool plusMinus) {
    // list of 1 and -1 representation

    if (plusMinus) {
        vector<int> S = genPlusMinus();
        assert(S.size() == A.size());

        vector<int> bestS = S;
        long bestResidue = plusMinusSum(S, A);
        long SResidue = bestResidue;

        for (int i = 0; i < g_NumIterations; i++) {
            vector<int> S_p = randMovePM(S, 5);
            long S_pResidue = plusMinusSum(S_p, A);
            long resDiff = -1 * (S_pResidue - SResidue);

            if(S_pResidue < SResidue) {
                S = S_p;
                SResidue = S_pResidue;
            } 
            // even occurs with probability prob_eT
            else if (rand()/RAND_MAX <= prob_eT(i, resDiff)) {
               S = S_p;
               SResidue = S_pResidue;
            }

            if (SResidue < bestResidue) {
                bestResidue = SResidue;
                bestS = S;
            }
        }
        return bestResidue;
    }
    // prepartition representation
    else {
        vector<int> P = genRandPartition();
        assert(P.size() == A.size());
        
        MaxHeap A_p = generateA_p(A, P);
        vector<int> bestP = P;
        long bestResidue = KarmarkarKarp(A_p);
        long prevResidue = bestResidue;

        for (int i = 0; i < g_NumIterations; i++) {
            vector<int> P_new = randMovePrepartition(P, 5);
            A_p = generateA_p(A, P_new);
            long currResidue = KarmarkarKarp(A_p);
            long resDiff = -1 * (currResidue - prevResidue);

            if (currResidue < prevResidue) {
                P = P_new;
                prevResidue = currResidue;
            }
            else if (rand()/RAND_MAX <= prob_eT(i, resDiff)) {
                P = P_new;
                prevResidue = currResidue;
            }

            if (currResidue < bestResidue) {
                bestResidue = currResidue;
                bestP = P_new;
            }
        }

        return bestResidue;
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        cout << "Usage: ./partition flag" << endl;
        return -1;
    }
    
    // 'A' is our list of values that we want to partition
    MaxHeap A;
    srand(time(NULL));

    // if flag is 0 we run the real algorithm
    if (strtol(argv[1], nullptr, 0) == 0) {
        if (argc != 4) {
            cout << "Usage: ./partition 0 algorithm inputfile" << endl;
            cout << "Algorithm options: 0, 1, 2, 3, 11, 12, 13" << endl;
            return -1;
        } 

        string filename(argv[3]);
        ifstream infile(filename);
        long entry;
        // insert all values from text file into heap A
        while (infile >> entry) {
            A.insert(entry);
        }

        N = A.size();

        int alg = strtol(argv[2], nullptr, 0);

        switch (alg) {
            case 0:
                cout << KarmarkarKarp(A) << endl;
                break;

            case 1:
                cout << repeatedRandom(A, true) << endl;
                break;

            case 2:
                cout << hillClimbing(A, true) << endl;
                break;

            case 3:
                cout << simulatedAnnealing(A, true) << endl;
                break;

            case 11:
                cout << repeatedRandom(A, false) << endl;
                break;
            
            case 12:
                cout << hillClimbing(A, false) << endl;
                break;

            case 13:
                cout << simulatedAnnealing(A, false) << endl;
                break;

            default:
                cout << "No algorithm for "  << alg << endl;
        }

        return 0;
    }

    // testing a single algorithm
    if (strtol(argv[1], nullptr, 0) == 1) {
        if (argc != 2) {
            cout << "Usage: ./partition 1" << endl;
            return -1;
        }

        vector<long> totalResidues(7, 0);
        vector<long> minResidue(7, LONG_MAX);

        for (int i = 0; i < 50; i++) {
            populateA(A);
            MaxHeap A_parent = A;

            long residue = KarmarkarKarp(A);
            totalResidues[0] += residue;
            A = A_parent;

            residue = repeatedRandom(A, true);
            totalResidues[1] += residue;
            if (residue < minResidue[1]) {
                minResidue[1] = residue;
            }
            A = A_parent;

            residue = hillClimbing(A, true);
            totalResidues[2] += residue;
            if (residue < minResidue[2]) {
                minResidue[2] = residue;
            }
            A = A_parent;

            residue = simulatedAnnealing(A, true);
            totalResidues[3] += residue;
            if (residue < minResidue[3]) {
                minResidue[3] = residue;
            }
            A = A_parent;

            residue = repeatedRandom(A, false);
            totalResidues[4] += residue;
            if (residue < minResidue[4]) {
                minResidue[4] = residue;
            }
            A = A_parent;

            residue = hillClimbing(A, false);
            totalResidues[5] += residue;
            if (residue < minResidue[5]) {
                minResidue[5] = residue;
            }
            A = A_parent;

            residue = simulatedAnnealing(A, false);
            totalResidues[6] += residue;
            if (residue < minResidue[6]) {
                minResidue[6] = residue;
            }
    
            A.clear();
        }

        cout << fixed;
        cout << setprecision(2);

        cout << "KK Algorithm: " << totalResidues[0]/50. << endl;

        cout << "Repeated Random PM: " << totalResidues[1]/50. << endl;
        cout << "min: " << minResidue[1] << endl;

        cout << "Hill Climbing PM: " << totalResidues[2]/50. << endl;
        cout << "min: " << minResidue[2] << endl;

        cout << "Simulated Annealing PM: " << totalResidues[3]/50. << endl;
        cout << "min: " << minResidue[3] << endl;

        cout << "Repeated Random Prepartitioned: " << totalResidues[4]/50. << endl;
        cout << "min: " << minResidue[4] << endl;

        cout << "Hill Climbing Prepartitioned: " << totalResidues[5]/50. << endl;
        cout << "min: " << minResidue[5] << endl;

        cout << "Simulated Annealing Prepartitioned: " << totalResidues[6]/50. << endl;
        cout << "min: " << minResidue[6] << endl;
    }

    return 0;
}



void MaxHeap::swap(int idx1, int idx2) {
    // swap locations in heap vector
    long buffer = heap[idx1];
    heap[idx1] = heap[idx2];
    heap[idx2] = buffer;
}

void MaxHeap::siftUp (int idx) {
    // if we are at root we're done
    if (idx == 0) {
        return;
    }
    int parIdx = (idx-1)/2;
    // if weight is greater than parent's weight bubble up
    if (heap[idx] > heap[parIdx]) {
        swap(idx, parIdx);
        siftUp(parIdx);
    }
}

void MaxHeap::siftDown (int idx) {
    int lChild = (2*idx) + 1;
    int rChild = (2*idx) + 2;

    int largest = idx;
    // find largest value index
    if (lChild < int(heapSize) && heap[lChild] > heap[idx]) {
        largest = lChild;
    }
    if (rChild < int(heapSize) && heap[rChild] > heap[largest]) {
        largest = rChild;
    }
    // if current idx is not smallest index, bubble down
    if (largest != idx)
    {
        swap(idx, largest);
        siftDown(largest);
    }

}

void MaxHeap::insert(long v) {
    heap.push_back(v);
    int idx = heapSize;
    heapSize++;
    siftUp(idx);
}

long MaxHeap::getMax () {
    // if heap is empty return nothing
    if (heapSize == 0) {
        // cout << "Heap is empty." << endl;
        return 0;
    }
    long toReturn = heap[0];
    // if we are removing the last element
    if (heapSize == 1) {
        heap.pop_back();
        heapSize--;
        return toReturn;
    }
    heap[0] = heap[heapSize-1];
    heap.pop_back();
    heapSize--;
    siftDown(0);
    return toReturn;
}

size_t MaxHeap::size() {
    return heap.size();
}

bool MaxHeap::empty() {
    return heapSize == 0;
}

void MaxHeap::clear() {
    heapSize = 0;
    heap.clear();
    heap.resize(0);
}