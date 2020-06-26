//
// Created by francesco on 31/05/20.
//

#ifndef GENETICTSP_SRC_GRAPH_GRAPH_H_
#define GENETICTSP_SRC_GRAPH_GRAPH_H_

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <iostream>
#include <random>

template<typename Key = int, typename Value = double>
class Graph {
 private:
  std::unordered_map<Key, std::unordered_map<Key, Value>> nodes;
 public:
  explicit Graph(int nodesNumber, const std::string &filepath = "none");
  Value GetEdgeValue(Key startNode, Key endNode);
  int GetNodesNumber() const;
  auto GetMapIterator() const;
};

/*** Construct the graph
 *
 * @tparam Key
 * @tparam Value
 * @param nodesNumber
 * @param filepath
 */
template<typename Key, typename Value>
Graph<Key, Value>::Graph(int nodesNumber, const std::string &filepath) {
  nodes.reserve(nodesNumber);

  if (filepath == "none") {
    std::random_device rd;    // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());   //Standard mersenne_twister_engine seeded with rd()
    std::uniform_real_distribution<double> unif(0, 100);
    //std::default_random_engine edgeValueGenerator;

    for (int index = 1; index <= nodesNumber; index++) {

      //! initialize first level of map
      nodes.insert(std::make_pair( index, std::unordered_map<Key, Value>()));
      nodes[index].reserve(nodesNumber - (index + 1));

      //! fill inner map
      for (int currentIndex = index + 1; currentIndex <= nodesNumber; currentIndex++) {
        nodes[index].insert(std::make_pair(currentIndex, unif(gen)));
      }
    }
    /*
    //! print map for test purpose
    for (auto exIt = nodes.begin(); exIt != nodes.end(); exIt++) {
      std::cout << "First key is: " << exIt->first << std::endl;
      for (auto innerIt = exIt->second.begin(); innerIt != exIt->second.end(); innerIt++) {
        std::cout << "Key: " << innerIt->first
                  << " Value:  " << innerIt->second << std::endl;
      }
      std::cout << std::endl;
    }*/

  } else {
    // parse file and construct the graph
  }
}

/*** Return the value of the edge connecting two node by
 *   First retrieve the outer map with [startNode], then
 *   take the edge value with [endNode]
 */
template<typename Key, typename Value>
Value Graph<Key, Value>::GetEdgeValue(Key startNode, Key endNode) {
  // case is a complete graph (check both nodes to obtain edge value)
  return (nodes[startNode].count(endNode)) ? nodes[startNode][endNode] : nodes[endNode][startNode];
}

template<typename Key, typename Value>
int Graph<Key, Value>::GetNodesNumber() const {
  return nodes.size();
}

/** Return a const iterator to iterate over
 *  map keys
 */
template<typename Key, typename Value>
auto Graph<Key, Value>::GetMapIterator() const {
  return nodes.cbegin();
}
#endif //GENETICTSP_SRC_GRAPH_GRAPH_H_