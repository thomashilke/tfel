#ifndef TIMER_H
#define TIMER_H

#include <chrono>


class timer {
public:
  timer(): start(std::chrono::high_resolution_clock::now()) {}
  double tic() {
    std::chrono::time_point<std::chrono::high_resolution_clock>
      now(std::chrono::high_resolution_clock::now());

    const std::chrono::duration<double> elapsed(now - start);
    return std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count();
  }

private:
  std::chrono::time_point<std::chrono::high_resolution_clock> start;
};

#endif /* TIMER_H */
