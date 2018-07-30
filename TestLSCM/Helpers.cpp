#include "Helpers.h"

#define NOMINMAX
#include <windows.h>

std::wstring s2ws (const std::string& s)
{
  int len;
  int slength = (int)s.length () + 1;
  len = MultiByteToWideChar (CP_ACP, 0, s.c_str (), slength, 0, 0);
  wchar_t* buf = new wchar_t[len];
  MultiByteToWideChar (CP_ACP, 0, s.c_str (), slength, buf, len);
  std::wstring r (buf);
  delete[] buf;
  return r;
}



//void DebugOut (const std::string fmt, ...)
//{
//  va_list args;
//  va_start (args, fmt);
//  size_t size = snprintf (nullptr, 0, fmt.c_str (), args) + 1;
//  // Extra space for '\0'
//  std::unique_ptr<char[]> buf (new char[size]);
//  snprintf (buf.get (), size, fmt.c_str (), args);
//  auto output = std::string (buf.get (), buf.get () + size - 1); // We don't want the '\0' inside
//
//  va_end (args);
//
//  output = output.append ("\r\n");
//
//  OutputDebugString ((LPCWSTR)s2ws (output).c_str ());
//}

