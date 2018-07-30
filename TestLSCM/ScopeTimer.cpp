#include "ScopeTimer.h"
#include "Helpers.h"



ScopeTimer::ScopeTimer ( const std::string& message )
  : _message(message)
  , _start( std::chrono::high_resolution_clock::now () )
{
}


ScopeTimer::~ScopeTimer ()
{
  auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now () - _start).count ();
  DebugOut ( _message + " in %lu milliseconds", ms );
}
