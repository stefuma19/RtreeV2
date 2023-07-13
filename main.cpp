#include <iostream>
#include "RTree.h"
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

typedef int ValueType;

#define DIM 2

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


int main() {
    //
    typedef RTree<ValueType, double, DIM, float, 4> MyTree;
    MyTree tree;

    //std::string filePath = "datasets/dataset_small.csv";
    std::string filePath = "../datasets/dataset_small.csv";

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

    /*
    auto list = tree.ListTree();
    int counter = 0;
    for (auto aabb: list) {
        cout << "TreeList [" << counter++ << "]: ";
        for (double j : aabb.m_min) {
            cout << j << " ";
        }
        cout << "; ";
        for (double j : aabb.m_max){
            cout << j << " ";
        }
        cout << endl;
    }
    */

    /*TODO: vedere implementazione di std::vector<typename RTREE_QUAL::Rect> RTREE_QUAL::ListTree() const per capire
    come iterare tra i vari nodi e avere accesso ai nodi interni */
    std::vector<double> query {0.5, 0.5};
    int k = 10;
    tree.myQuery(k, query);

    return 0;

}