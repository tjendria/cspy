#ifndef SRC_CC_SEARCH_H__
#define SRC_CC_SEARCH_H__

#include <memory> // unique_ptr
#include <set>
#include <vector>

#include "labelling.h"

namespace bidirectional {

class Search {
 public:
  /// Direction of Search
  Directions direction;
  //  Params     params;
  /// Stopping criteria for each direction
  bool stop = false;
  /// Stopping criteria for each direction
  bool bound_exceeded = false;
  /// Number of unprocessed labels generated
  int unprocessed_count = 0;
  /// Number of labels processed
  int processed_count = 0;
  /// Number of labels generated (includes the possibly infeasible extensions)
  int generated_count = 0;

  size_t k;

  /* Search-related parameters */

  /// Lower bounds from any node to sink
  std::unique_ptr<std::vector<double>> lower_bound_weight;
  /// vector with indices of vertices visited
  std::set<int>                     visited_vertices;
  std::shared_ptr<labelling::Label> current_label;
  /// Intermediate current best label with possibly complete source-sink path
  /// (shared pointer as we want to be able to substitute it without
  /// resetting)

  // std::shared_ptr<labelling::Label> intermediate_label;
  /// Need custom comparator for multiset since elements are pointers
  struct LabelCompare {
    using CT = std::shared_ptr<labelling::Label>;
    bool operator()(const CT& a, const CT& b) const {
      const int c_res = a->params_ptr->critical_res;
      auto getScalar = [c_res](const CT& label) {
        return (label->resource_consumption.size() > c_res && // push invalid labels to end
                label->partial_path.size() > 1)
                        ? label->resource_consumption[c_res]
                        : std::numeric_limits<double>::max();
      };
      return getScalar(a) < getScalar(b);
    }
  };
  /// Store k-best intermediate labels, use multiset to keep them ordered
  std::multiset<std::shared_ptr<labelling::Label>, LabelCompare>
      intermediate_labels;

  /// vector with pareto optimal labels (per node) in each direction
  std::vector<std::vector<labelling::Label>> efficient_labels;
  /// vector with pointer to label with least weight (per node) in each
  /// direction
  std::vector<std::shared_ptr<labelling::Label>> best_labels;
  /**
   * heap vector to keep unprocessed labels ordered.
   * the order depends on the on the direction of the search.
   * i.e. forward -> increasing in the monotone resource,
   * backward -> decreasing in the monotone resource.
   */
  std::unique_ptr<std::vector<labelling::Label>> unprocessed_labels;

  // TODO: Use bucket-heap
  /* Heap operations for vector of labels */

  /**
   * Initialises heap using the appropriate comparison
   * i.e. increasing in the monotone resource forward lists, decreasing
   * otherwise
   */
  void makeHeap();

  /**
   * Push new elements in heap using the appropriate comparison
   * i.e. increasing in the monotone resource forward lists, decreasing
   * otherwise
   */
  void pushHeap();

  void pushUnprocessedLabel(const labelling::Label& label) {
    unprocessed_labels->push_back(label);
    pushHeap();
  }

  void pushEfficientLabel(const int& lemon_id, const labelling::Label& label) {
    efficient_labels[lemon_id].push_back(label);
  }

  /// Replace best label
  void replaceBestLabel(const int& lemon_id, const labelling::Label& label) {
    auto label_ptr = std::make_shared<labelling::Label>(label);
    best_labels[lemon_id].swap(label_ptr);
  }

  /// Replace best label
  void replaceCurrentLabel(const labelling::Label& label) {
    auto label_ptr = std::make_shared<labelling::Label>(label);
    current_label.swap(label_ptr);
  }

  /// Replace intermediate label
  void replaceIntermediateLabel(const labelling::Label& label) {
    auto label_ptr = std::make_shared<labelling::Label>(label);

    // Remove the current intermediate label if it is the only one
    if (intermediate_labels.size() == 1 && (*intermediate_labels.begin())->partial_path.size() == 1) {
      intermediate_labels.clear();
    }

    intermediate_labels.insert(label_ptr);

    while (intermediate_labels.size() >= k) {
      intermediate_labels.erase(*intermediate_labels.rbegin());
    }
  }

  /// Update vertices visited
  void addVisitedVertex(const int& lemon_id) {
    visited_vertices.insert(lemon_id);
  }

  Search(const Directions& direction_in, size_t k = 4);
  ~Search(){};
};

} // namespace bidirectional

#endif // SRC_CC_SEARCH_H__
