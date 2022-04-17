#include <iostream>
#include <iomanip>
#include <vector>
// #include <utility>
#include <fstream>
#include <time.h>
#include <stdlib.h>
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
        long val2 = A.getMax();
        // if val2 is 0 then all other entries are 0, so 
        //      val1 is the residue
        if (val2 == 0) {
            return val1;
        }

        long difference = val1 - val2;

        A.insert(difference);
        A.insert(0);
    }
}

// calculate residue when partition is represented as list S of 1 and -1
long plusMinusSum(vector<int> S, MaxHeap &A) {
    long residue = 0;
    assert(S.size() == A.size());

    for (int i=0; i < S.size(); i++) {
        residue += A.heap[i]*S[i];
    }
    return abs(residue);
}

vector<int> randMovePM(vector<int> &S) {
    vector<int> S_p = S;
    srand(time(NULL));
    int i = rand() % N;
    // set si to -si with probability 1/2
    if (rand() % 2 == 0) {
        S_p[i] = -1 * S_p[i];
    }
    // generate j such that j != i
    int j = rand() % N;
    while (i == j) {
        j = rand() % N;
    }
    // set sj = -sj
    S_p[j] = -1 * S_p[j];

    return S_p;
}

vector<int> randMovePrepartition(vector<int> &P) {
    vector<int> P_new = P;
    srand(time(NULL));
    int i = rand() % N;
    // generate j such that j != i
    int j = rand() % N;
    while (P_new[i] == j) {
        j = rand() % N;
    }
    // set p_i = j
    P_new[i] = j;

    return P_new;
}

// given a prepartition, generate A' using algorithm as described in pset
MaxHeap generateA_p(MaxHeap &A, vector<int> &P) {
    vector<long> V_p(N, 0);
    MaxHeap A_p;

    for (int j = 0; j < N; j++) {
        int p_j = P[j];
        V_p[p_j] += V_p[p_j] + A.heap[j];
    }

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
    srand(time(NULL));
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
    srand(time(NULL));
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
            if (residue > bestResidue) {
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
            long residue = KarmarkarKarp(A_p);
            if (residue > bestResidue) {
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
            vector<int> S_p = randMovePM(S);
            long currResidue = plusMinusSum(S_p, A);
            if (currResidue < bestResidue) {
                bestResidue = currResidue;
                bestS = S_p;
            }
            S = S_p;
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
            vector<int> P_new = randMovePrepartition(P);
            A_p = generateA_p(A, P_new);
            long currResidue = KarmarkarKarp(A_p);

            if (currResidue < bestResidue) {
                bestResidue = currResidue;
                bestP = P_new;
            }
            P = P_new;
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
            vector<int> S_p = randMovePM(S);
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
        
        vector<int> bestP = P;
        MaxHeap A_p = generateA_p(A, P);
        long bestResidue = KarmarkarKarp(A_p);

        for (int i = 0; i < g_NumIterations; i++) {
            vector<int> P_new = randMovePrepartition(P);
            A_p = generateA_p(A, P_new);
            long currResidue = KarmarkarKarp(A_p);

            if (currResidue < bestResidue) {
                bestResidue = currResidue;
                bestP = P_new;
            }
            P = P_new;
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
    }

    // testing a single algorithm
    if (strtol(argv[1], nullptr, 0) == 1) {
        if (argc != 2) {
            cout << "Usage: ./partition 1" << endl;
            return -1;
        }

        vector<long> totalResidues(7, 0);

        for (int i = 0; i < 50; i++) {
            populateA(A);
            MaxHeap A_parent = A;

            long residue = KarmarkarKarp(A);
            totalResidues[0] += residue;
            A = A_parent;

            residue = repeatedRandom(A, true);
            totalResidues[1] += residue;
            A = A_parent;

            residue = hillClimbing(A, true);
            totalResidues[2] += residue;
            A = A_parent;

            residue = simulatedAnnealing(A, true);
            totalResidues[3] += residue;
            A = A_parent;

            residue = repeatedRandom(A, false);
            totalResidues[4] += residue;
            A = A_parent;

            residue = hillClimbing(A, false);
            totalResidues[5] += residue;
            A = A_parent;

            residue = simulatedAnnealing(A, false);
            totalResidues[6] += residue;
    
            A.clear();
        }

        cout << fixed;
        cout << setprecision(2);

        cout << "KK Algorithm: " << totalResidues[0]/50. << endl;
        cout << "Repeated Random PM: " << totalResidues[1]/50. << endl;
        cout << "Hill Climbing PM: " << totalResidues[2]/50. << endl;
        cout << "Simulated Annealing PM: " << totalResidues[3]/50. << endl;
        cout << "Repeated Random Prepartitioned: " << totalResidues[4]/50. << endl;
        cout << "Hill Climbing Prepartitioned: " << totalResidues[5]/50. << endl;
        cout << "Simulated Annealing Prepartitioned: " << totalResidues[6]/50. << endl;
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
    if (lChild < heapSize && heap[lChild] > heap[idx]) {
        largest = lChild;
    }
    if (rChild < heapSize && heap[rChild] > heap[largest]) {
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
    siftUp(idx);
    heapSize++;
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