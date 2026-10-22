#pragma once

#include <cstddef>
#include <exception>
#include <vector>

namespace zSpace { namespace zTsDetail {

// One slot per iteration, allocated before entering OpenMP. Never let an
// exception escape a worker's structured block; rethrow after the join.
// The lowest failing index wins, independently of the scheduling order.
class ParallelErrors {
  std::vector<std::exception_ptr> errors;
public:
  explicit ParallelErrors(std::size_t count) : errors(count) {}
  void capture(std::size_t index) noexcept {
    errors[index] = std::current_exception();
  }
  void rethrow() const {
    for (const auto &error : errors)
      if (error) std::rethrow_exception(error);
  }
};

}} // namespace zSpace::zTsDetail
