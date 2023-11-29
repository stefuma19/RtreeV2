#include "RTree.h"
#include <iostream>
#include <utility>
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <nlopt.hpp>

using namespace std;

typedef int ValueType;
typedef std::vector<double> MyTuple;

//Helper class to replace decimal separator
template<class charT, charT sep>

class punct_facet : public std::numpunct<charT> {
protected:
    charT do_decimal_point() const override { return sep; }
};


#define DIM 2      //Dimensions of the dataset
#define ONLY_RTREE //If defined, only the r-tree method is executed, the sequential ones are ignored
#define BETA 0.66  //Beta value for the directional query

#if DIM == 2       //Define the leaf capacity depending on the number of dimensions
#define LEAF_CAPACITY 100
#endif
#if DIM == 3
#define LEAF_CAPACITY 70
#endif
#if DIM == 4
#define LEAF_CAPACITY 55
#endif
#if DIM == 5
#define LEAF_CAPACITY 45
#endif
#if DIM == 6
#define LEAF_CAPACITY 39
#endif

struct Rect {
    Rect() = default;

    // Constructor for n dimensions
    Rect(const std::vector<double> &a_min, const std::vector<double> &a_max) {
        min = a_min;
        max = a_max;
    }

    std::vector<double> min;
    std::vector<double> max;
};

// Define the Object class
struct NodeWithScore {

    explicit NodeWithScore(double s) {
        score = s;
    }

    double score;

    bool operator<(const NodeWithScore &a) const {
        return score < a.score;
    }
};


bool compareLastColumn(const MyTuple &tuple1, const MyTuple &tuple2);

bool compareLastColumn(const MyTuple &tuple1, const MyTuple &tuple2) {
    // Assuming the last column is accessed using tuple.back()
    return tuple1.back() < tuple2.back();
}

std::vector<Rect> createRectanglesFromCSV(const string &filename, int numDimensions) {
    std::vector<Rect> rectangles;
    std::ifstream file(filename);

    if (!file) {
        std::cout << "Failed to open the file: " << filename << std::endl;
        return rectangles;
    }

    std::string line;
    std::getline(file, line); // Skip the header line

    while (std::getline(file, line)) {
        std::vector<double> minValues(numDimensions);
        std::vector<double> maxValues(numDimensions);

        std::stringstream ss(line);
        std::string value;
        int dimension = 0;

        // Skip the first column
        std::getline(ss, value, ',');

        while (std::getline(ss, value, ',')) {
            //cout << score << endl;
            if (dimension < numDimensions) {
                minValues[dimension] = std::stod(value);
                maxValues[dimension] = std::stod(value);
            }
            dimension++;
        }

        if (dimension == numDimensions) {
            rectangles.emplace_back(minValues, maxValues);
        } else {
            std::cout << "Invalid row format: " << line << std::endl;
        }
    }

    file.close();

    return rectangles;
}

double *convertDoubleVectorToIntPointer(const vector<double> &doubleVec) {
    static std::vector<double> intVec; // Static vector to ensure the memory is not deallocated

    // Convert std::vector<double> to std::vector<int>
    intVec.assign(doubleVec.begin(), doubleVec.end());

    // Obtain a pointer to the first element of the std::vector<int>
    double *intPtr = intVec.data();

    return intPtr;
}

// Compute the distance between a point and the preference line given the query vector. The preference line is computed
// as the inverse of the query vector inside this function.
double dist_line_point(const vector<double> &point, const vector<double> &prefLine, double den) {

    double dist = 0.0;

    for (int i = 0; i < DIM; i++) {
        double num = 0.0;

        for (int j = 0; j < DIM; j++) {
            num += prefLine[j] * point[j];
        }

        double d = point[i] - (prefLine[i] * num / den);
        d = d * d;
        dist += d;
    }

    return std::sqrt(dist);
}

//Creates a data vector containing the tuples of the dataset from the given csv file
std::vector<MyTuple> createDataVectorFromCSV(const string &filename) {

    std::vector<MyTuple> data;  // Vector to store the vectors of values

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
    }
    std::string line;
    std::getline(file, line);  // Skip the first line (header)

    while (std::getline(file, line)) {
        std::vector<double> values; //Current tuple
        std::istringstream iss(line);
        std::string token;
        // Skip the first value in each line (the id)
        std::getline(iss, token, ',');

        // Read and save the remaining values in values
        while (std::getline(iss, token, ',')) {
            values.push_back(std::stod(token));
        }
        data.push_back(values); // Add current tuple to vector of tuples
    }
    file.close();

    return data;
}

//Performs a linear top-k query given the dataset saved in a vector
void sequentialLinearWithVector(std::vector<MyTuple> points, std::vector<double> query, int k) {

    std::priority_queue<NodeWithScore> heap;

    double score;
    int size = points.size();
    int i;
    int j;

    auto startTimeLinSeq = std::chrono::high_resolution_clock::now();

    for (i = 0; i < k; i++) {
        score = 0;
        for (j = 0; j < DIM; j++) {
            score += points[i][j] * query[j]; //Compute the score of the current point
        }
        heap.push(NodeWithScore(score)); //Insert the point in the heap
    }

    for (; i < size; i++) {
        score = 0;
        for (j = 0; j < DIM; j++) {
            score += points[i][j] * query[j]; //Compute the score of the current point
        }

        if (score < heap.top().score) {
            heap.pop();
            heap.emplace(score); //Insert the point in the heap if its score is greater than the worst score in the heap
        }
    }

    auto endTimeLinSeq = std::chrono::high_resolution_clock::now();
    auto durationLinSeq = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeLinSeq - startTimeLinSeq);
    std::cout << "[K = " << k << " ] - Linear Sequential: " << durationLinSeq.count() << " ms" << std::endl;
}

// Reads the CSV file containing the dataset and performs the sequential linear top-k query
void sequentialLinear(const std::string &filename, std::vector<double> query, int k) {

    std::priority_queue<NodeWithScore> heap;

    std::vector<MyTuple> data;  // Vector to store the vectors of values
    double score;
    int i;
    int j = 0;
    //Current line
    auto startTimeLinSeq = std::chrono::high_resolution_clock::now();
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
    }
    std::string line;
    std::getline(file, line);  // Skip the first line (header)
    while (std::getline(file, line)) {
        std::vector<double> values;
        std::istringstream iss(line);
        std::string token;
        score = 0;
        // Skip the first value in each line
        std::getline(iss, token, ',');

        i = 0;
        // Read the remaining values
        while (std::getline(iss, token, ',')) {
            score += std::stod(token) * query[i];
            i++;
            values.push_back(std::stod(token));
            i++;
        }

        data.push_back(values);

        if (j < k) {
            heap.push(NodeWithScore(score));
        } else if (score < heap.top().score) {
            heap.pop();
            heap.push(NodeWithScore(score));
        }
        j++;
    }
    file.close();
    auto endTimeLinSeq = std::chrono::high_resolution_clock::now();
    auto durationLinSeq = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeLinSeq - startTimeLinSeq);
    std::cout << "Linear Sequential: " << durationLinSeq.count() << " ms" << std::endl;

    sequentialLinearWithVector(data, std::move(query), k);
}

//Performs sequential directional query score computation if data is already stored in a vector
void sequentialDirectionalWithVector(std::vector<MyTuple> points, std::vector<double> query, int k) {
    std::priority_queue<NodeWithScore> heap;

    double score;
    int size = points.size();
    int i;
    int j = 0;

    auto startTimeDirSeq = std::chrono::high_resolution_clock::now();

    std::vector<double> prefLine = computePreferenceLine(query);

    double den = 0;
    den = computeDenFromPrefLine(prefLine);

    for (i = 0; i < k; i++) { //The first k elements will certainly add in the initial top-k
        score = 0;
        for (j = 0; j < DIM; j++) {
            score += points[i][j] * query[j];
            score = BETA * score + (1 - BETA) * dist_line_point(points[i], prefLine, den);
        }
        heap.emplace(score);
    }

    for (; i < size; i++) {
        score = 0;
        for (j = 0; j < DIM; j++) {
            score += points[i][j] * query[j];
            score = BETA * score + (1 - BETA) * dist_line_point(points[i], prefLine, den);
        }

        if (score < heap.top().score) {
            heap.pop();
            heap.emplace(score);
        }
    }
    auto endTimeDirSeq = std::chrono::high_resolution_clock::now();
    auto durationDirSeq = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeDirSeq - startTimeDirSeq);
    std::cout << "[K = " << k << " ] - Directional Sequential: " << durationDirSeq.count() << " ms" << std::endl;

#ifdef PRINT_RESULTS
    for (i = 0; i < k; i++) {
        std::cout << heap.top().score << std::endl;
        heap.pop();
    }
#endif


}

/*
// Reads the CSV file containing the dataset and performs the sequential directional top-k query
void sequentialDirectional(const std::string &filename, std::vector<double> query, int k) {
    std::priority_queue<NodeWithScore> heap;

    std::vector<MyTuple> data;  // Vector to store the vectors of values
    double score;
    int i;
    int j = 0;

    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
    }
    auto startTimeDirSeq = std::chrono::high_resolution_clock::now();
    std::string line;
    std::getline(file, line);  // Skip the first line (header)

    double point[DIM];

    for (j = 0; j < k; j++) {
        std::getline(file, line);
        std::istringstream iss(line);
        std::string token;
        score = 0;
        // Skip the first value in each line
        std::getline(iss, token, ',');
        i = 0;
        // Read the remaining values
        while (std::getline(iss, token, ',')) {
            point[i] = (std::stod(token));
            score += std::stod(token) * query[i];
            i++;
        }

        score = BETA * score + (1 - BETA) * dist_line_point(point, query);
        heap.push(NodeWithScore(score));
    }

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        score = 0;
        // Skip the first value in each line
        std::getline(iss, token, ',');
        i = 0;
        // Read the remaining values
        while (std::getline(iss, token, ',')) {
            point[i] = (std::stod(token));
            score += std::stod(token) * query[i];
            i++;
        }

        score = BETA * score + (1 - BETA) * dist_line_point(point, query);

        if (score < heap.top().score) {
            heap.pop();
            heap.emplace(score);
        }
    }

    file.close();
    auto endTimeDirSeq = std::chrono::high_resolution_clock::now();
    auto durationDirSeq = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeDirSeq - startTimeDirSeq);
    std::cout << "Directional Sequential: " << durationDirSeq.count() << " ms" << std::endl;

}
*/

int main() {

    typedef RTree<ValueType, double, DIM, float, LEAF_CAPACITY> MyTree;
    MyTree tree;

    // PATH CONTAINING THE DATASET TO PROCESS
    std::string datasetPath = "../datasets/cor_neg/2D/cor_neg_1M_2.csv";


    #ifndef ONLY_RTREE
    std::vector<double> queryIniziale = {0.2, 0.8};
    std::vector<int> k_values = {1, 5, 10, 50, 100, 500};

    std::vector<MyTuple> data;  // Vector to store the vectors of values
    data = createDataVectorFromCSV(datasetPath);

    for (int k_value: k_values) {
        //printf("K = %d\n", k_value);
        //sequentialLinear(datasetPath, queryIniziale, k_value);
        //sequentialDirectional(datasetPath, queryIniziale, k_value);
        sequentialLinearWithVector(data, queryIniziale, k_value);
        sequentialDirectionalWithVector(data, queryIniziale, k_value);
    }

    //return 0;
    /*
    sequentialLinear(datasetPath, queryIniziale, 10);
    sequentialDirectional(datasetPath, queryIniziale, 10);
    */
    #endif

    //Tree creation
    std::vector<Rect> rectangles = createRectanglesFromCSV(datasetPath, DIM);

    int i = 0;
    // Process the created rectangles as needed
    for (const Rect &rect: rectangles) {
        std::vector<double> minValues = rect.min;
        std::vector<double> maxValues = rect.max;

        double *min = convertDoubleVectorToIntPointer(minValues);
        double *max = convertDoubleVectorToIntPointer(minValues);

        tree.Insert(min, max, i);
        i++;
    }

    std::cout << "----------------R_TREE STATS----------------" << std::endl;
    int numNodes, numLeaves, numPoints;
    tree.CountBoxesLeavesAndPoints(&numNodes, &numLeaves, &numPoints);
    std::cout << "contBox: " << numNodes << std::endl;
    std::cout << "contLeaf: " << numLeaves << std::endl;
    std::cout << "contPoint: " << numPoints << std::endl;

    int k = 10;
    int numQ = 0;

    double timeLinRT = 0;
    double timeDirRT = 0;
    double timeLinSeq = 0;
    double timeDirSeq = 0;

    int numBoxLin = 0;
    int numLeavesLin = 0;
    int numPointLin = 0;

    int numBoxDir = 0;
    int numLeavesDir = 0;
    int numPointDir = 0;

    std::vector<double> timeLinRTvect;
    timeLinRTvect.reserve(9);
    std::vector<double> timeDirRTvect;
    timeLinRTvect.reserve(9);

    std::vector<double> numBoxLinvect;
    timeLinRTvect.reserve(9);
    std::vector<double> numLeavesLinvect;
    timeLinRTvect.reserve(9);
    std::vector<double> numPointLinvect;
    timeLinRTvect.reserve(9);

    std::vector<double> numBoxDirvect;
    timeLinRTvect.reserve(9);
    std::vector<double> numLeavesDirvect;
    timeLinRTvect.reserve(9);
    std::vector<double> numPointDirvect;
    timeLinRTvect.reserve(9);

    std::vector<int> k_values;

    std::vector<double> totalNonLinearProblemsSolvedVect;
    std::vector<double> totalNonLinearProblemsSolvedTimeVect;

    std::vector<double> totalLinearScoresComputedVect;
    std::vector<double> totalTimeLinearScoresComputationVect;

    std::vector<double> totalDirectionalScoresComputedVect;
    std::vector<double> totalTimeDirectionalScoresComputationVect;

    std::vector<double> totalLinearBoundsComputedVect;
    std::vector<double> totalTimeLinearBoundsComputationVect;

    std::vector<double> totalLinearScoresComputedForDirectionalVect;
    std::vector<double> totalTimeLinearScoresComputationForDirectionalVect;

    std::vector<double> totalLinearBoundsComputedForDirectionalVect;
    std::vector<double> totalTimeLinearBoundsComputationForDirectionalVect;


    std::ifstream kValuesFile("../utilities/new_k.txt"); //FILE CONTAINING THE K VALUES TO PROCESS

    std::cout << "\n - Results for dataset = " << datasetPath << std::endl;

    while (kValuesFile >> k) {

        //std::cout << "\n - Results for k = " << k << ", dataset = " << datasetPath << std::endl;

        timeLinRT = 0;
        timeDirRT = 0;
        timeLinSeq = 0;
        timeDirSeq = 0;

        numBoxLin = 0;
        numLeavesLin = 0;
        numPointLin = 0;

        numBoxDir = 0;
        numLeavesDir = 0;
        numPointDir = 0;

        std::ifstream queriesFile("../queries/" + std::to_string(DIM) + "d.txt"); //FILE CONTAINING THE QUERIES TO PROCESS
        //std::ifstream queriesFile("../balanced_queries/" + std::to_string(DIM) + "d.txt");
        //std::ifstream queriesFile("../unbalanced_queries/" + std::to_string(DIM) + "d.txt");
        numQ = 0;

        if (queriesFile.is_open()) {
            std::string line;
            while (std::getline(queriesFile, line)) { // Read each line of the file
                numQ++;
                std::vector<double> query;
                std::stringstream ss(line);
                std::string value;

                while (std::getline(ss, value, ',')) { // Split line by comma
                    double number = std::stod(value);
                    query.push_back(number);
                }

                //std::cout << "----------------LINEAR RTREE----------------" << std::endl;
                auto startTimeLinRT = std::chrono::high_resolution_clock::now();
                tree.linearTopKQueryRTree(k, query, &numBoxLin, &numLeavesLin, &numPointLin);
                auto endTimeLinRT = std::chrono::high_resolution_clock::now();
                auto durationLinRT = std::chrono::duration_cast<std::chrono::microseconds>(
                        endTimeLinRT - startTimeLinRT);

                timeLinRT += durationLinRT.count();

                //std::cout << "----------------DIRECTIONAL RTREE----------------" << std::endl;
                auto startTimeDirRT = std::chrono::high_resolution_clock::now();

                #if TREE_METHOD == TIGHT
                tree.DirectionalTopKQueryRTree(k, query, &numBoxDir, &numLeavesDir, &numPointDir);
                #endif

                #if TREE_METHOD == LOOSE
                tree.DirectionalTopKQueryRTreeLoose(k, query, &numBoxDir, &numLeavesDir, &numPointDir);
                #endif
                //tree.DirectionalTopKQueryRTreeMixed(k, query, &numBoxDir, &numLeavesDir, &numPointDir);
                auto endTimeDirRT = std::chrono::high_resolution_clock::now();
                auto durationDirRT = std::chrono::duration_cast<std::chrono::microseconds>(
                        endTimeDirRT - startTimeDirRT);

                timeDirRT += durationDirRT.count();
            }
        }

        //Print to have a useful output for the Excel
        timeLinRTvect.push_back(static_cast<double> (timeLinRT) / numQ);
        timeDirRTvect.push_back(static_cast<double> (timeDirRT) / numQ);

        numBoxLinvect.push_back(static_cast<double>(numBoxLin) / numQ);
        numLeavesLinvect.push_back(static_cast<double>(numLeavesLin) / numQ);
        numPointLinvect.push_back(static_cast<double>(numPointLin) / numQ);

        numBoxDirvect.push_back(static_cast<double>(numBoxDir) / numQ);
        numLeavesDirvect.push_back(static_cast<double>(numLeavesDir) / numQ);
        numPointDirvect.push_back(static_cast<double>(numPointDir) / numQ);

        #ifdef MEASURE_TIME

        k_values.push_back(k);

        /* NON LINEAR PROBLEMS */
        totalNonLinearProblemsSolvedVect.push_back(static_cast<double>(totalNonLinearProblemsSolved)/numQ);
        totalNonLinearProblemsSolvedTimeVect.push_back(static_cast<double>(totalTimeNonLinearProblemsExecution) /numQ);
        //RESET totalNonLinearProblemsSolved and time
        totalNonLinearProblemsSolved = 0;
        totalTimeNonLinearProblemsExecution = 0;

        /* LINEAR SCORES (LINEAR ONLY)*/
        totalLinearScoresComputedVect.push_back(static_cast<double>(totalLinearScoresComputed)/numQ);
        totalTimeLinearScoresComputationVect.push_back(static_cast<double>(totalTimeLinearScoresComputation)/numQ);

        totalLinearScoresComputed = 0;
        totalTimeLinearScoresComputation = 0;

        /* DIRECTIONAL SCORES */
        totalDirectionalScoresComputedVect.push_back(static_cast<double>(totalDirectionalScoresComputed)/numQ);
        totalTimeDirectionalScoresComputationVect.push_back(static_cast<double>(totalTimeDirectionalScoresComputation)/numQ);

        totalDirectionalScoresComputed = 0;
        totalTimeDirectionalScoresComputation = 0;

        /* LINEAR BOUNDS (LINEAR ONLY) */
        totalLinearBoundsComputedVect.push_back(static_cast<double>(totalLinearBoundsComputed)/numQ);
        totalTimeLinearBoundsComputationVect.push_back(static_cast<double>(totalTimeLinearBoundsComputation)/numQ);

        totalLinearBoundsComputed = 0;
        totalTimeLinearBoundsComputation = 0;

        /* LINEAR SCORES (FOR DIRECTIONAL) */
        totalLinearScoresComputedForDirectionalVect.push_back(static_cast<double>(totalLinearScoresComputedForDirectional)/numQ);
        totalTimeLinearScoresComputationForDirectionalVect.push_back(static_cast<double>(totalTimeLinearScoresComputationForDirectional)/numQ);

        totalLinearScoresComputedForDirectional = 0;
        totalTimeLinearScoresComputationForDirectional = 0;

        /* LINEAR BOUNDS (FOR DIRECTIONAL) */
        totalLinearBoundsComputedForDirectionalVect.push_back(static_cast<double>(totalLinearBoundsComputedForDirectional)/numQ);
        totalTimeLinearBoundsComputationForDirectionalVect.push_back(static_cast<double>(totalTimeLinearBoundsComputationForDirectional)/numQ);

        totalLinearBoundsComputedForDirectional = 0;
        totalTimeLinearBoundsComputationForDirectional = 0;

        #endif
    }

    //Print to have a useful output for the Excel
    std::cout << "\n - Results for Rtree Execution" << std::endl;

    /*
    std::cout << "\nLinear Rtree time [us]:" << std::endl;
    for (const auto &element: timeLinRTvect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }

    std::cout << "\nLinear Rtree numPoint:" << std::endl;
    for (const auto &element: numPointLinvect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }
    std::cout << "\nLinear Rtree numLeaves:" << std::endl;
    for (const auto &element: numLeavesLinvect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }
    std::cout << "\nLinear Rtree numBoxes:" << std::endl;
    for (const auto &element: numBoxLinvect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }

    #ifdef MEASURE_TIME
    std::cout << "\nLinear Rtree Bounds Computed:" << std::endl;
    for (const auto &element: totalLinearBoundsComputedVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }
    std::cout << "\nLinear Rtree Bounds Computation Time: [us]" << std::endl;
    for (const auto &element: totalTimeLinearBoundsComputationVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element/1000 << std::endl; //raw data is in ns and needs to be printed in us
    }

    std::cout << "\nLinear Rtree Scores:" << std::endl;
    for (const auto &element: totalLinearScoresComputedVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }
    std::cout << "\nLinear Rtree Scores Computation Time: [us]" << std::endl;
    for (const auto &element: totalTimeLinearScoresComputationVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element/1000 << std::endl; //raw data is in ns and needs to be printed in us
    }

    #endif

    std::cout << "\nDirectional Rtree time [us]:" << std::endl;
    for (const auto &element: timeDirRTvect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }

    std::cout << "\nDirectional Rtree numPoint:" << std::endl;
    for (const auto &element: numPointDirvect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }
    std::cout << "\nDirectional Rtree numLeaves:" << std::endl;
    for (const auto &element: numLeavesDirvect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }
    std::cout << "\nDirectional Rtree numBoxes:" << std::endl;
    for (const auto &element: numBoxDirvect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }

    #ifdef MEASURE_TIME
    std::cout << "\nDirectional Rtree Problems:" << std::endl;
    for (const auto &element: totalNonLinearProblemsSolvedVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }
    std::cout << "\nDirectional Rtree ProblemTime [us]:" << std::endl;
    for (const auto &element: totalNonLinearProblemsSolvedTimeVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element/1000 << std::endl; // raw data is in ns and needs to be printed in us
    }

    std::cout << "\nDirectional Rtree Scores:" << std::endl;
    for (const auto &element: totalDirectionalScoresComputedVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }
    std::cout << "\nDirectional Rtree Scores Computation Time: [us]" << std::endl;
    for (const auto &element: totalTimeDirectionalScoresComputationVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element/1000 << std::endl; //raw data is in ns and needs to be printed in us
    }

    std::cout << "\nLinear Rtree Scores Computed For Directional:" << std::endl;
    for (const auto &element: totalLinearScoresComputedForDirectionalVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }
    std::cout << "\nLinear Rtree Scores Computed For Directional Computation Time: [us]" << std::endl;
    for (const auto &element: totalTimeLinearScoresComputationForDirectionalVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element/1000 << std::endl; //raw data is in ns and needs to be printed in us
    }

    std::cout << "\nLinear Rtree Bounds Computed For Directional:" << std::endl;
    for (const auto &element: totalLinearBoundsComputedForDirectionalVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element << std::endl;
    }
    std::cout << "\nLinear Rtree Bounds Computed For Directional Computation Time: [us]" << std::endl;
    for (const auto &element: totalTimeLinearBoundsComputationForDirectionalVect) {
        std::cout.imbue(std::locale(std::cout.getloc(), new punct_facet<char, ','>));
        std::cout << element/1000 << std::endl; //raw data is in ns and needs to be printed in us
    }
    #endif

    */

    #ifdef MEASURE_TIME

    // Store vectors in a vector of vectors as strings
    std::vector<std::vector<std::string>> allVectors;
    std::vector<std::string> firstRow;

    firstRow.push_back("n");
    firstRow.push_back("d");
    firstRow.push_back("k");
    firstRow.push_back("Rtree boxes");
    firstRow.push_back("Rtree leaves");
    firstRow.push_back("Method");
    firstRow.push_back("Linear Rtree time [us]");
    firstRow.push_back("Linear Rtree numPoint");
    firstRow.push_back("Linear Rtree numLeaves");
    firstRow.push_back("Linear Rtree numBoxes");
    firstRow.push_back("Linear Rtree Bounds Computed");
    firstRow.push_back("Linear Rtree Bounds Computation Time: [us]");
    firstRow.push_back("Linear Rtree Scores");
    firstRow.push_back("Linear Rtree Scores Computation Time: [us]");
    firstRow.push_back("Directional Rtree time [us]");
    firstRow.push_back("Directional Rtree numPoint");
    firstRow.push_back("Directional Rtree numLeaves");
    firstRow.push_back("Directional Rtree numBoxes");
    firstRow.push_back("Directional Rtree Problems");
    firstRow.push_back("Directional Rtree ProblemTime [us]");
    firstRow.push_back("Directional Rtree Scores");
    firstRow.push_back("Directional Rtree Scores Computation Time: [us]");
    //firstRow.push_back("Linear Rtree Scores Computed For Directional");
    //firstRow.push_back("Linear Rtree Scores Computed For Directional Computation Time: [us]");
    firstRow.push_back("Linear Rtree Bounds Computed For Directional");
    firstRow.push_back("Linear Rtree Bounds Computed For Directional Computation Time: [us]");
    allVectors.push_back(firstRow);
    for (size_t i = 0; i < totalTimeLinearBoundsComputationForDirectionalVect.size(); ++i) {
        std::vector<std::string> row;
        row.push_back(std::to_string(numPoints));
        row.push_back(std::to_string(DIM));
        row.push_back(std::to_string(k_values[i]));
        row.push_back(std::to_string(numNodes));
        row.push_back(std::to_string(numLeaves));
        #if TREE_METHOD == TIGHT
        row.push_back("Tight");
        #endif
        #if TREE_METHOD == LOOSE
        row.push_back("Loose");
        #endif
        row.push_back(std::to_string(timeLinRTvect[i]));
        row.push_back(std::to_string(numPointLinvect[i]));
        row.push_back(std::to_string(numLeavesLinvect[i]));
        row.push_back(std::to_string(numBoxLinvect[i]));
        row.push_back(std::to_string(totalLinearBoundsComputedVect[i]));
        row.push_back(std::to_string(totalTimeLinearBoundsComputationVect[i]/1000));
        row.push_back(std::to_string(totalLinearScoresComputedVect[i]));
        row.push_back(std::to_string(totalTimeLinearScoresComputationVect[i]/1000));
        row.push_back(std::to_string(timeDirRTvect[i]));
        row.push_back(std::to_string(numPointDirvect[i]));
        row.push_back(std::to_string(numLeavesDirvect[i]));
        row.push_back(std::to_string(numBoxDirvect[i]));
        row.push_back(std::to_string(totalNonLinearProblemsSolvedVect[i]));
        row.push_back(std::to_string(totalNonLinearProblemsSolvedTimeVect[i]/1000));
        row.push_back(std::to_string(totalDirectionalScoresComputedVect[i]));
        row.push_back(std::to_string(totalTimeDirectionalScoresComputationVect[i]/1000));
        //row.push_back(std::to_string(totalLinearScoresComputedForDirectionalVect[i]));
        //row.push_back(std::to_string(totalTimeLinearScoresComputationForDirectionalVect[i]/1000));
        row.push_back(std::to_string(totalLinearBoundsComputedForDirectionalVect[i]));
        row.push_back(std::to_string(totalTimeLinearBoundsComputationForDirectionalVect[i]/1000));
        allVectors.push_back(row);
    }

    // Printing vectors as CSV to the console
    /*
    for (const auto& row : allVectors) {
        for (size_t i = 0; i < row.size(); ++i) {
            std::cout << row[i];

            if (i != row.size() - 1) {
                std::cout << ",";
            }
        }
        std::cout << "\n";
    }
     */

    // Open a file for writing
    std::ofstream outFile("../stats.csv");

    // Check if the file is open
    if (outFile.is_open()) {

        for (const auto& row : allVectors) {
            for (size_t t = 0; t < row.size(); ++t) {
                outFile << row[t];

                if (t != row.size() - 1) {
                    outFile << ",";
                }
            }
            outFile << "\n";
        }

        // Close the file
        outFile.close();
        std::cout << "File written successfully." << std::endl;
    } else {
        std::cerr << "Unable to open the file." << std::endl;
    }

    #endif
    return 0;
}