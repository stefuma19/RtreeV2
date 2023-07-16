#include <iostream>
#include "RTree.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <chrono>
#include <nlopt.hpp>

using namespace std;

typedef int ValueType;
typedef std::vector<double> MyTuple;


#define DIM 2
#define BETA 0.66


struct Rect
{
    Rect() {}

    // Constructor for n dimensions
    Rect(const std::vector<double>& a_min, const std::vector<double>& a_max)
    {
        min = a_min;
        max = a_max;
    }

    std::vector<double> min;
    std::vector<double> max;
};

bool compareLastColumn(const MyTuple& tuple1, const MyTuple& tuple2) {
    // Assuming the last column is accessed using tuple.back()
    return tuple1.back() < tuple2.back();
}

std::vector<Rect> createRectanglesFromCSV(const std::string& filename, int numDimensions)
{
    std::vector<Rect> rectangles;
    std::ifstream file(filename);

    if (!file)
    {
        std::cout << "Failed to open the file: " << filename << std::endl;
        return rectangles;
    }

    std::string line;
    std::getline(file, line); // Skip the header line

    while (std::getline(file, line))
    {
        std::vector<double> minValues(numDimensions);
        std::vector<double> maxValues(numDimensions);

        std::stringstream ss(line);
        std::string value;
        int dimension = 0;

        // Skip the first column
        std::getline(ss, value, ',');

        while (std::getline(ss, value, ','))
        {
            //cout << score << endl;
            if (dimension < numDimensions){
                minValues[dimension] = std::stod(value);
                maxValues[dimension] = std::stod(value);
            }
            dimension++;
        }

        if (dimension == numDimensions)
        {
            rectangles.emplace_back(minValues, maxValues);
        }
        else
        {
            std::cout << "Invalid row format: " << line << std::endl;
        }
    }

    file.close();

    return rectangles;
}

double* convertDoubleVectorToIntPointer(const std::vector<double>& doubleVec)
{
    static std::vector<double> intVec; // Static vector to ensure the memory is not deallocated

    // Convert std::vector<double> to std::vector<int>
    intVec.assign(doubleVec.begin(), doubleVec.end());

    // Obtain a pointer to the first element of the std::vector<int>
    double* intPtr = intVec.data();

    return intPtr;
}

double dist_line_point(const std::vector<double>& point, const std::vector<double>& query) {
    double dist = 0.0;
    int len = query.size();

    for (int i = 0; i < len; i++) {
        double num = 0.0;
        double den = 0.0;

        for (int j = 0; j < len; j++) {
            num += query[j] * point[j];
            den += query[j] * query[j];
        }

        double d = point[i] - (query[i] * num / den);
        d = d * d;
        dist += d;
    }

    return std::sqrt(dist);
}

std::vector<MyTuple> readCSVLin(const std::string& filename, std::vector<double> query) {
    std::vector<MyTuple> data;  // Vector to store the vectors of values
    double score;
    int i;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return data;  // Return an empty vector
    }
    std::string line;
    std::getline(file, line);  // Skip the first line (header)
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        std::vector<double> values;
        score = 0;
        // Skip the first value in each line
        std::getline(iss, token, ',');

        i= 0;
        // Read the remaining values
        while (std::getline(iss, token, ',')) {
            values.push_back(std::stod(token));
            score += std::stod(token) * query[i];
            i++;
        }
        values.push_back(score);
        // Add the vector of values to the vector
        data.push_back(values);
    }

    file.close();
    return data;
}

std::vector<MyTuple> readCSVDir(const std::string& filename, std::vector<double> query) {
    std::vector<MyTuple> data;  // Vector to store the vectors of values
    double score;
    int i;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return data;  // Return an empty vector
    }
    std::string line;
    std::getline(file, line);  // Skip the first line (header)
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string token;
        std::vector<double> values;
        score = 0;
        // Skip the first value in each line
        std::getline(iss, token, ',');

        i= 0;
        // Read the remaining values
        while (std::getline(iss, token, ',')) {
            values.push_back(std::stod(token));
            score += std::stod(token) * query[i];
            i++;
        }

        std::vector<double> point = {values.begin(), values.end()};
        score = BETA * score + (1 - BETA) * dist_line_point(point, query);
        values.push_back(score);
        // Add the vector of values to the vector
        data.push_back(values);
    }

    file.close();
    return data;
}


int main() {
    //
    typedef RTree<ValueType, double, DIM, float, 100> MyTree;
    MyTree tree;

    //std::string filePath = "../datasets/dataset_small.csv";
    std::string filePath = "../datasets/cor_neg_1k_2.csv";
    //std::string filePath = "../datasets/cor_neg_1M_2.csv";
    //std::string filePath = "../datasets/cor_neg_1M_4.csv";
    //std::string filePath = "../datasets/cor_neg_1k_4.csv";

    //Tree creation
    std::vector<Rect> rectangles = createRectanglesFromCSV(filePath, DIM);

    int i = 0;
    // Process the created rectangles as needed
    for (const Rect& rect : rectangles)
    {
        std::vector<double> minValues = rect.min;
        std::vector<double> maxValues = rect.max;

        double* min = convertDoubleVectorToIntPointer(minValues);
        double* max = convertDoubleVectorToIntPointer(minValues);

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

    std::ifstream inputFile("../utilities/k.txt");

    while (inputFile >> k) {

        std::cout << "\n - Results for k = " << k << ", dataset = " << filePath << std::endl;

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

        std::ifstream file("../queries/" + std::to_string(DIM) + "d.txt");
        numQ = 0;

        if (file.is_open()) {
            std::string line;
            while (std::getline(file, line)) { // Read each line of the file
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
                tree.DirectionalTopKQueryRTree(k, query, &numBoxDir, &numLeavesDir, &numPointDir);
                auto endTimeDirRT = std::chrono::high_resolution_clock::now();
                auto durationDirRT = std::chrono::duration_cast<std::chrono::microseconds>(
                        endTimeDirRT - startTimeDirRT);

                timeDirRT += durationDirRT.count();
            }
        }

        std::cout << "\nExecution time LinRT: " << static_cast<double> (timeLinRT) / numQ << " microseconds." << std::endl;
        std::cout << "Execution time DirRT: " << static_cast<double>(timeDirRT) / numQ << " microseconds." << std::endl;

        std::cout << "\nAccesses to Box Linear: " << static_cast<double>(numBoxLin) / numQ << std::endl;
        std::cout << "Accesses to Leaves Linear: " << static_cast<double>(numLeavesLin) / numQ << std::endl;
        std::cout << "Accesses to Points Linear: " << static_cast<double>(numPointLin) / numQ << std::endl;

        std::cout << "\nAccesses to Box Directional: " << static_cast<double>(numBoxDir) / numQ << std::endl;
        std::cout << "Accesses to Leaves Directional: " << static_cast<double>(numLeavesDir) / numQ << std::endl;
        std::cout << "Accesses to Points Directional: " << static_cast<double>(numPointDir) / numQ << std::endl;
    }

    std::vector<double> query;
    query.reserve(DIM);
    for (int i = 0; i < DIM; i++) {
        query.push_back(1 / DIM);
    }

    //std::cout << "----------------LINEAR SEQUENTIAL----------------" << std::endl;
    auto startTimeLinSeq = std::chrono::high_resolution_clock::now();

    std::vector<MyTuple> tuplesLin = readCSVLin(filePath, query);
    std::sort(tuplesLin.begin(), tuplesLin.end(), compareLastColumn);

    auto endTimeLinSeq = std::chrono::high_resolution_clock::now();
    auto durationLinSeq = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeLinSeq - startTimeLinSeq);


    /*for (int i = k-1; i >= 0; i--){
        std::cout << i << " tuples score: " << tuplesLin[i].back() << std::endl;
    }*/

    timeLinSeq += durationLinSeq.count();

    //std::cout << "----------------DIRECTIONAL SEQUENTIAL----------------" << std::endl;
    auto startTimeDirSeq = std::chrono::high_resolution_clock::now();

    std::vector<MyTuple> tuplesDir = readCSVDir(filePath, query);
    std::sort(tuplesDir.begin(), tuplesDir.end(), compareLastColumn);

    auto endTimeDirSeq = std::chrono::high_resolution_clock::now();
    auto durationDirSeq = std::chrono::duration_cast<std::chrono::milliseconds>(endTimeDirSeq - startTimeDirSeq);

    /*for (i = k-1; i >= 0; i--){
        std::cout << i << " tuples score: " << tuplesDir[i].back() << std::endl;
    }*/

    timeDirSeq += durationDirSeq.count();

    std::cout << "\nResults for Sequential Execution" << std::endl;
    std::cout << "\nExecution time LinSeq: " << timeLinSeq << " milliseconds." << std::endl;
    std::cout << "Execution time DirSeq: " << timeDirSeq << " milliseconds." << std::endl;

    return 0;
}