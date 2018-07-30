#pragma once

#include <memory>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdio>


#define NOMINMAX
#include <Windows.h>

#define DEL( ptr ) \
if (ptr)              \
{                     \
    delete ptr;       \
    ptr = nullptr;    \
}

#define DEL_ARR( ptr ) \
if (ptr)              \
{                     \
    delete [] ptr;    \
    ptr = nullptr;    \
}

std::wstring s2ws (const std::string& s);

template<typename ... Args>
void DebugOut (const std::string fmt, Args ... args)
{
#if _DEBUG  
  auto output = Format (fmt, args ...);

  output = output.append ("\r\n");

  //OutputDebugString ((LPCWSTR)s2ws (output).c_str ());

  std::cout << output;
#endif
}

template<typename ... Args>
std::string Format (const std::string fmt, Args ... args)
{
  size_t size = snprintf (nullptr, 0, fmt.c_str (), args ...) + 1;
  // Extra space for '\0'
  std::unique_ptr<char[]> buf (new char[size]);
  snprintf (buf.get (), size, fmt.c_str (), args ...);
  return std::string (buf.get (), buf.get () + size - 1); // We don't want the '\0' inside
}







