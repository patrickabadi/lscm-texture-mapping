#pragma once

#include <string>
#include <chrono>

class ScopeTimer
{
public:
  ScopeTimer (const std::string& message);
  ~ScopeTimer ();

private:
  std::string _message;
  std::chrono::time_point<std::chrono::steady_clock> _start;
};

